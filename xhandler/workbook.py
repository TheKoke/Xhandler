import os
import numpy as np

from analysis import *


class WorkbookParser:
    def __init__(self, path: str, spectres_path: str) -> None:
        self.path = path
        self.spectres_path = spectres_path

        self.all = None

    def find_angle(self, angle: float) -> Spectrum:
        return next(analyzed for analyzed in self.all if analyzed.angle == angle)

    def collect_all(self) -> list[Spectrum]:
        if 'workbook.txt' not in os.listdir(self.path):
            return []

        all_reports = open(self.path + '\\workbook.txt', 'r').read().split('\n\n')[1:]
        return [self.reverse_parameters(report) for report in all_reports]

    def reverse_parameters(self, report: str) -> Spectrum:
        result = Spectrum()

        result.angle = self.__get_angle()
        result.scale_value, result.scale_shift = self.__get_calibration_constants(report)
        result.peaks = self.__get_peaks(report)

        return result

    def __get_spectrum(self, path: str, angle: float) -> np.ndarray:
        return np.loadtxt(path + f'{int(angle)}.txt')

    def __get_angle(self, report: str) -> float:
        fiducial_index = report.index('angle') - 2
        return float(report[:fiducial_index])

    def __get_calibration_constants(self, report: str) -> tuple[float, float]:
        equation = self.__get_equation(report)
        parts = equation.split(' + ')

        scale_shift = float(parts[1])
        scale_value = float(parts[0].split(' * ')[0])

        return (scale_value, scale_shift)

    def __get_equation(self, report: str) -> str:
        lines = report.split('\n')
        info_str = 'Calibrated by equation: E(ch) = '

        position = lines[1].index(info_str) + len(info_str)
        return lines[1][position + 1:]

    #TODO: Implement this method.
    def __get_peaks(self, report: str) -> list[PeakSupervisor]:
        lines = report.split('\n')
        info_str = 'Peaks analysis info'

        peaks_start = next(i for i in range(len(lines)) if info_str in lines[i]) + 2
        taken = []

        for i in range(peaks_start, len(lines)):
            pass

        return taken
    
    # TODO: Implement this method.
    def __reverse_peak(self, center: float, fwhm: float, calib_coeff: float, calib_e0: float) -> PeakSupervisor:
        fwtm = fwhm / np.log10(2)

        channel_width = int(fwtm / calib_coeff)
        center_in_channel = int((center - calib_e0) / calib_coeff)

        spectrum = self.__get_spectrum()


# TODO: Croutch.
class WorkbookMaster:
    def __init__(self, spectres_path: str) -> None:
        self.path = os.getcwd()
        self.spectres_path = spectres_path

        self.parser = WorkbookParser(self.path, self.spectres_path)

    def write(self, message: str) -> None:
        file = open(self.path + '\\workbook.txt', 'a')
        file.write(message)
        file.close()

    def gather_analyzed(self) -> list[Analytics]:
        return self.parser.all


if __name__ == '__main__':
    pass
