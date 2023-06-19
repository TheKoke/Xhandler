import numpy as np
from global_const import MASS_EXCESSES, STATES


class Nuclei:
    def __init__(self, nuclons: int, charge: int) -> None:
        self.nuclons = nuclons
        self.charge = charge

    def __str__(self) -> str:
        return f'A: {self.nuclons}, Z: {self.charge}'

    @property
    def mass_excess(self) -> float:
        return MASS_EXCESSES[(self.charge, self.nuclons)]
    
    @property
    def states(self) -> list[float]:
        return STATES[(self.charge, self.nuclons)]
    
    @property
    def wigner_width(self) -> list[float]:
        pass
    
    def mass(self, unit: str = 'MeV') -> float:
        if unit == 'MeV':
            return self.charge * 938.27 + (self.nuclons - self.charge) * 939.57
        
        if unit == 'a.m.u.':
            return self.charge * 1.007276467 + (self.nuclons - self.charge) * 1.008664915
        
        if unit == 'g':
            return self.charge * 1.672e-24 + (self.nuclons - self.charge) * 1.675e-24


class Reaction:
    def __init__(self, beam: Nuclei, target: Nuclei, fragment: Nuclei, beam_energy: float) -> None:
        self.beam = beam
        self.target = target
        self.fragment = fragment
        self.residual = self.__residual_nuclei()

        self.beam_energy = beam_energy

    @property
    def is_elastic(self) -> bool:
        return self.beam == self.fragment

    def __residual_nuclei(self) -> Nuclei:
        nuclon = (self.beam.nuclons + self.target.nuclons) - self.fragment.nuclons
        charge = (self.beam.charge + self.target.charge) - self.fragment.charge
        return Nuclei(nuclon, charge)

    def reaction_quit(self, residual_state: float = 0) -> float:
        q0 = (self.beam.mass_excess + self.target.mass_excess) - (self.fragment.mass_excess + self.residual.mass_excess)
        return q0 - residual_state
    
    def fragment_energy(self, residual_state: float, fragment_angle: float) -> np.ndarray:
        r = Reaction.__r_factor(
            self.beam.mass(), 
            self.beam_energy, 
            self.fragment.mass(), 
            self.residual.mass(), 
            fragment_angle * np.pi / 180
        )

        s = Reaction.__s_factor(
            self.beam.mass(), 
            self.beam_energy, 
            self.fragment.mass(), 
            self.residual.mass(), 
            self.reaction_quit(residual_state)
        )

        return (r + np.sqrt(r ** 2 + s)) ** 2
    
    def residual_energy(self, residual_state: float) -> np.ndarray:
        r = Reaction.__r_factor(
            self.beam.mass(),
            self.beam_energy, 
            self.residual.mass(),
            self.fragment.mass(),
            self.residual_angle(residual_state)
        )

        s = Reaction.__s_factor(
            self.beam.mass(),
            self.beam_energy,
            self.residual.mass(),
            self.fragment.mass(),
            self.reaction_quit(residual_state)
        )

        return (r + np.sqrt(r ** 2 + s)) ** 2
    
    def residual_angle(self, residual_state: float, fragment_angle: float) -> float:
        fragment_ears = self.fragment_energy(residual_state)
        energy_relation = np.sqrt(self.beam.mass() * self.beam_energy / (self.fragment.mass() * fragment_ears))

        return np.pi / 2 - np.arctan(
            (energy_relation - np.cos(fragment_angle * np.pi / 180)) / np.sin(fragment_angle * np.pi / 180)
        )
    
    @staticmethod
    def __r_factor(beam_mass: float, beam_energy: float, 
                   instance_mass: float, partner_mass: float, angle: float) -> float:
        numerator = np.sqrt(beam_mass * instance_mass * beam_energy) * np.cos(angle)
        return numerator / (instance_mass + partner_mass)

    @staticmethod
    def __s_factor(beam_mass: float, beam_energy: float, 
                   instance_mass: float, partner_mass: float, reaction_quit: float) -> float:
        numerator = beam_energy * (partner_mass - beam_mass) + partner_mass * reaction_quit
        return numerator / (instance_mass + partner_mass)


class Ionization:
    def __init__(self, stray: Nuclei, environ: Nuclei, detector: str = 'Ge') -> None:
        self.stray = stray
        self.environ = environ
        self.detector = detector

    def energy_loss(self, energy: float, thickness: float) -> None:
        return self.specific_energy_loss(energy) * thickness
    
    def specific_energy_loss(self, energy: float) -> float:
        electron_mass = 0.511 # MeV
        reduced_planck = 6.582e-22 # MeV * s
        lightspeed = 3e10 # cm / s
        fine_structure = 1 / 137 # dimensionless

        e_power_4 = (reduced_planck * lightspeed * fine_structure) ** 2
        betta_power_2 = self.lorenz_parameter(energy) ** 2

        common = 4 * np.pi * self.electrons_density() * self.stray.charge ** 2 * e_power_4 / (electron_mass * betta_power_2)
        logarithm = np.log(2 * electron_mass * betta_power_2 / self.mean_environ_excitation())
        relativistic = np.log(1 - betta_power_2) + betta_power_2

        return common * (logarithm - relativistic)

    def mean_environ_excitation(self) -> float:
        hydrogen_ionization = 13.6e-6 # MeV
        return hydrogen_ionization * self.environ.charge
    
    def electrons_density(self) -> float:
        ro = self.__set_ro() # g * sm^-3

        avogadro = 6.02e23 # mol^-1
        return self.environ.charge * ro * avogadro / self.environ.nuclons # electrons * sm^-3
    
    def __set_ro(self) -> float:
        if self.detector.lower() == 'ge':
            return 5.323 # g * sm^-3
        if self.detector.lower() == 'si':
            return 2.33 # g * sm^-3
        
        return 337.8e-6
    
    def lorenz_parameter(self, energy: float) -> float:
        #  dimensionless     MeV          MeV        
        return np.sqrt(2 * energy / self.stray.mass())


if __name__ == '__main__':
    pass
