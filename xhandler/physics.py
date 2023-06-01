import numpy as np


MASS_EXCESSES = {
    (1, 1)  : 7.28900, (1, 2)  : 13.1357, (1, 3)  : 14.9498,
    (2, 3)  : 14.9312, (2, 4)  : 2.42490,
    (3, 6)  : 14.0869, (3, 7)  : 14.9071, (3, 8)  : 20.9458,
    (4, 7)  : 15.7690, (4, 9)  : 11.3485, (4, 10) : 12.6075, (4, 11) : 20.1772,
    (5, 8)  : 22.9216, (5, 10) : 12.0506, (5, 11) : 8.66770, (5, 12) : 13.3694,
    (6, 10) : 15.6987, (6, 11) : 10.6494, (6, 12) : 0.00000, (6, 13) : 3.12500, (6, 14) : 3.01990,
    (7, 12) : 17.3381, (7, 13) : 5.34550, (7, 14) : 2.86340, (7, 15) : 0.10140, (7, 16) : 5.68390,
    (8, 14) : 8.00780, (8, 15) : 2.85560, (8, 16) : -4.7370, (8, 17) : -0.8088, (8, 18) : -0.7828, (8, 19) : 3.3329,
    (9, 17) : 1.95170, (9, 18) : 0.87310, (9, 19) : -1.4874, (9, 20) : -0.0175, (9, 21) : -0.0476, (9, 22) : 2.7934,
}


class Nuclei:
    def __init__(self, nuclons: int, charge: int) -> None:
        self.nuclons = nuclons
        self.charge = charge
        self.mass_excess = self.__mass_excess()

    def __str__(self) -> str:
        return f'A: {self.nuclons}, Z: {self.charge}'

    def __mass_excess(self) -> float:
        return MASS_EXCESSES[(self.charge, self.nuclons)]

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
    

if __name__ == '__main__':
    pass
