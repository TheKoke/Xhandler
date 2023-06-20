import numpy as np


class Peak:
    def __init__(self) -> None:
        self.__mu = 0
        self.__fwhm = 1
        self.__area = 1

    @property
    def mu(self) -> float:
        return self.__mu
    
    @mu.setter
    def mu(self, input: float) -> None:
        self.__mu = input

    @property
    def fwhm(self) -> float:
        return self.__fwhm
    
    @fwhm.setter
    def fwhm(self, input: float) -> None:
        if input <= 0:
            raise ValueError("FWHM of peak should be greater than zero.")
        self.__fwhm = input

    @property
    def area(self) -> float:
        return self.__area
    
    @area.setter
    def area(self, input: float) -> None:
        if input <= 0:
            raise ValueError("Area of peak should be greater than zero.")
        self.__area = input

    @property
    def hwhm(self) -> float:
        return self.__fwhm / 2


class Spectrum:
    def __init__(self) -> None:
        self.__angle = 0
        self.__spectrum = np.array([])

        self.__scale_value = 0
        self.__scale_shift = 0

        self.__peaks: list[Peak] = list()

    @property
    def angle(self) -> float:
        return self.__angle
    
    @angle.setter
    def angle(self, pretend: float) -> None:
        if pretend <= 0 or pretend > 180:
            raise ValueError("Reaction angle must be in interval: 0 < angle <= 180")
        self.__angle = pretend

    @property
    def spectrum(self) -> np.ndarray:
        return self.__spectrum
    
    @spectrum.setter
    def spectrum(self, input: np.ndarray) -> None:
        if len(input) == 0:
            raise ValueError("Spectrum doesn't be empty")
        self.__spectrum = input

    @property
    def energy_view(self) -> np.ndarray:
        return np.arange(1, len(self.__spectrum) + 1) * self.scale_value + self.scale_shift

    @property
    def calibration(self) -> tuple[float, float]:
        return (self.scale_value, self.scale_shift)

    @property
    def scale_value(self) -> float:
        return self.__scale_value

    @scale_value.setter
    def scale_value(self, input: float) -> None:
        if input <= 0:
            raise ValueError("Scale value can't be negative or equal to 0.")
        self.__scale_value = input

    @property
    def scale_shift(self) -> float:
        return self.__scale_shift
    
    @scale_shift.setter
    def scale_shift(self, input: float) -> None:
        self.__scale_shift = input

    @property
    def is_calibrated(self) -> bool:
        return not self.__scale_value <= 0 and not self.__scale_shift == 0

    @property
    def peaks(self) -> list[Peak]:
        return self.__peaks
    
    @peaks.setter
    def peaks(self, pretend: list[Peak]) -> None:
        self.__peaks = pretend


if __name__ == '__main__':
    pass
