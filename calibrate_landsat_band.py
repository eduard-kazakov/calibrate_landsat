class CalibrateLandsatBand():
    def __init__(self, band_file, metadata_file, band_number='0'):
        self.metadata_reader = LandsatMetadataReader(metadata_file)

        self.metadata = self.metadata_reader.metadata

        if band_number == '0':
            self.band_metadata = self.metadata_reader.get_band_metadata_by_file_name(band_file)
        else:
            self.band_metadata = self.metadata_reader.get_band_metadata_by_number(band_number)

        if not self.band_metadata:
            raise KeyError('Invalid band')


        self.band_dataset = gdal.Open(band_file)
        self.band_array = self.band_dataset.GetRasterBand(1).ReadAsArray()


    def get_radiance_as_array(self):
        radiance = ((self.band_metadata['radiance_maximum']-self.band_metadata['radiance_minimum']) / (self.band_metadata['quantize_cal_maximum']-self.band_metadata['quantize_cal_minimum'])) * (self.band_array - self.band_metadata['quantize_cal_minimum']) + self.band_metadata['radiance_minimum']
        radiance[self.band_array==0] = np.nan
        return radiance

    def get_reflectance_as_array(self, not_native_radiance_array=False):
        if self.band_metadata['type'] != 'reflectance':
            raise TypeError('Given band is thermal')
        if type(not_native_radiance_array)==bool:
            radiance = self.get_radiance_as_array()
        else:
            radiance = not_native_radiance_array
        d = float(self.metadata['EARTH_SUN_DISTANCE'])
        O = np.deg2rad(float(self.metadata['SUN_ELEVATION']))
        E = self.band_metadata['solar_irradiance']

        reflectance = (np.pi*radiance*d*d)/(E*np.sin(O))
        return reflectance

    def get_brightness_temperature_as_array(self):
        if self.band_metadata['type'] != 'thermal':
            raise TypeError('Given band is reflectance')

        radiance = self.get_radiance_as_array()

        K1 = self.band_metadata['k1_constant']
        K2 = self.band_metadata['k2_constant']

        brightness_temperature = (K2 / (np.log((K1/radiance+1)))) - 273.15

        return brightness_temperature

    def perform_dos_correction_for_radiance (self, radiance_array):
        dark_object_radiance = np.nanpercentile(radiance_array,0.01)
        O = np.deg2rad(float(self.metadata['SUN_ELEVATION']))
        E = self.band_metadata['solar_irradiance']
        d = float(self.metadata['EARTH_SUN_DISTANCE'])
        L_1p = (0.01*np.cos(O)*np.cos(O)*np.cos(O)*E) / (np.pi*d*d)
        Lhaze = dark_object_radiance - L_1p
        corrected_radiances = radiance_array - Lhaze
        corrected_reflectance = self.get_reflectance_as_array(not_native_radiance_array=corrected_radiances)
        corrected_reflectance[corrected_reflectance < 0] = 0
        return corrected_reflectance

    def calculate_physical_temperature_with_rad_transfer_model (self, lse_array, upWelling_radiance, downWelling_radiance, atmospheric_transmittance, ndvi_array=None):
        if self.band_metadata['type'] != 'thermal':
            raise TypeError('Given band is reflectance')

        # calc toaRadiance
        radiance = self.get_radiance_as_array()

        if ndvi_array != None:
            # calc LSE zhang angoritm
            condition_list = [np.logical_and(ndvi_array < -0.185, ndvi_array >= -1),
                             np.logical_and(ndvi_array >= -0.185, ndvi_array <= 0.157),
                             np.logical_and(ndvi_array >= 0.157, ndvi_array <= 0.727),
                             np.logical_and(ndvi_array > 0.727, ndvi_array <= 1)]

            mixed_pixels = np.log(ndvi_array)
            mixed_pixels = np.multiply(mixed_pixels, 0.047)
            mixed_pixels = np.add(mixed_pixels, 1.009)
            choice_list = [0.995, 0.985, mixed_pixels, 0.990]
            lse_array = np.select(condition_list, choice_list)

        # calc LST radiativeTransferEquation
        K1 = self.band_metadata['k1_constant']
        K2 = self.band_metadata['k2_constant']

        left = np.subtract(radiance, float(upWelling_radiance))
        te = np.multiply(float(atmospheric_transmittance), lse_array)
        left = np.divide(left, te)

        right = np.subtract(1, lse_array)

        cond_list = [lse_array < 0, lse_array > 0]
        choice_list = [1, np.divide(right, lse_array)]
        right = np.select(cond_list, choice_list)

        right = np.multiply(right, float(downWelling_radiance))
        LTS = np.subtract(left, right)

        plank_lower = np.divide(K1, LTS)
        plank_lower = np.add(plank_lower, 1)

        log_cond_list = [plank_lower == 0, lse_array > 0]
        log_choice_list = [1, np.log(plank_lower)]
        plank_lower = np.select(log_cond_list, log_choice_list)

        plank_cond_list = [plank_lower == 0, plank_lower > 0]
        plank_choice_list = ['nan', (np.divide(K2, plank_lower))]
        lst = np.select(plank_cond_list, plank_choice_list)
        return lst

    def save_array_as_gtiff(self, array, new_file_path):
        driver = gdal.GetDriverByName("GTiff")
        dataType = gdal.GDT_Float32
        dataset = driver.Create(new_file_path, self.band_dataset.RasterXSize, self.band_dataset.RasterYSize, self.band_dataset.RasterCount, dataType)
        dataset.SetProjection(self.band_dataset.GetProjection())
        dataset.SetGeoTransform(self.band_dataset.GetGeoTransform())
        dataset.GetRasterBand(1).WriteArray(array)
        del dataset
