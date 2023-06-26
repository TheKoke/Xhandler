import numpy as np


class Gaussian:
    def __init__(self, mu: float, fwhm: float, area: float) -> None:
        self.mu = mu
        self.fwhm = fwhm
        self.area = area

    def __str__(self) -> str:
        func = 'G(x) = '

        func += f'{np.round(self.area, 3)} / (sqrt(2pi * {np.round(self.dispersion(), 3)}) * '
        func += f'exp[-(x - {np.round(self.mu, 3)})^2 / 2{np.round(self.dispersion(), 3)}^2]'

        return func
    
    @property
    def dispersion(self) -> float:
        return self.fwhm / (2 * np.sqrt(2 * np.log(2)))

    def three_sigma(self) -> np.ndarray:
        return np.linspace(-3 * self.dispersion + self.mu, 3 * self.dispersion + self.mu, 50)

    def function(self) -> np.ndarray:
        constant = self.area / np.sqrt(2 * np.pi * self.dispersion ** 2)
        exp_constant = -1 / (2 * self.dispersion ** 2)

        array_part = (self.three_sigma() - self.mu) ** 2

        return constant * np.exp(exp_constant * array_part)
    

class Lorentzian:
    def __init__(self, mu: float, fwhm: float, area: float) -> None:
        self.mu = mu
        self.fwhm = fwhm
        self.area = area

    def __str__(self) -> str:
        func = 'L(x) = '
        func += f'2 * {round(self.area, 3)} / pi'
        func += f'* [{np.round(self.fwhm, 3)} / ((x - {round(self.mu)})^2 + {round(self.fwhm, 3)}^2)]'

        return func
    
    def dispersion(self) -> float:
        return self.fwhm / (2 * np.sqrt(2 * np.log(2)))
    
    def three_sigma(self) -> np.ndarray:
        sigma = self.dispersion()
        return np.linspace(-5 * sigma + self.mu, 5 * sigma + self.mu, 50)
    
    def pdf(self) -> np.ndarray:
        constant = 2 * self.area / np.pi
        brackets = self.fwhm / (4 * (self.three_sigma() - self.mu) ** 2 + self.fwhm ** 2)

        return constant * brackets


if __name__ == '__main__':
    pass
