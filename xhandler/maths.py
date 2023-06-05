import numpy as np


class Gaussian:
    def __init__(self, xy: np.ndarray) -> None:
        self.xdata = xy[0]
        self.ydata = xy[1]

        self.peak_index = len(self.xdata) // 2
        self.peak_center = self.xdata[self.peak_index]

        self.area, self.dispersion = self.__approximate()

    def __str__(self, /, dispersion_view: bool = False, fwhm_view: bool = True) -> str:
        func = 'G(x) = '

        if dispersion_view:
            func += f'{np.round(self.area, 3)} / (sqrt(2pi * {np.round(self.dispersion(), 3)}) * ' \
            f'exp(-(x - {np.round(self.peak_center, 3)})^2 / 2{np.round(self.dispersion(), 3)})'

        if fwhm_view:
            func += f'{np.round(self.area, 3)} / ({np.round(self.fwhm(), 3)} * sqrt(pi / 4ln2)) * ' \
            f'exp(- 4ln2 * (x - {np.round(self.peak_center, 3)})^2 / {np.round(self.fwhm() ** 2, 3)}))'

        return func
    
    def __approximate(self) -> tuple[float, float]:
        ls_params = self.__least_squares()

        dispersion = np.sqrt(-1 / (2 * ls_params[1]))
        area = np.sqrt(2 * np.pi * dispersion ** 2) * np.exp(ls_params[0])

        return (area, dispersion)

    def __least_squares(self) -> tuple[float, float]:
        xs = (self.xdata - self.peak_center) ** 2
        ys = np.log(self.ydata)

        coeff_matrix = np.array([[len(xs), xs.sum()], [xs.sum(), (xs ** 2).sum()]])
        right_side = np.array([ys.sum(), (xs * ys).sum()])

        return np.linalg.solve(coeff_matrix, right_side)

    def fwhm(self) -> np.float64:
        return self.dispersion * (2 * np.sqrt(2 * np.log(2)))

    def three_sigma(self) -> np.ndarray:
        return np.linspace(-3 * self.dispersion + self.peak_center, 3 * self.dispersion + self.peak_center, 50)

    def function(self) -> np.ndarray:
        constant = self.area / np.sqrt(2 * np.pi * self.dispersion ** 2)
        exp_constant = -1 / (self.dispersion ** 2)

        array_part = (self.three_sigma() - self.peak_center) ** 2

        return constant * np.exp(exp_constant * array_part)


if __name__ == '__main__':
    pass
