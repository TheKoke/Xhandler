from enum import Enum
from physics import Reaction, Nuclei
from ensdf import NAME2CHARGE


class ReactionNotation(Enum):
    SOVETIAN = 1
    CHEMIST = 2
    UNDEFINED = 3


class ReactionMaster:
    def __init__(self, input: str, energy: float) -> None:
        '''
        Nuclear reactions can wroted in 2 different styles:
        A(B, C)D - sovetian variant.
        A + B -> C + D - chemistry variant.
        '''
        self.reaction = input
        self.energy = energy

        self.notation = self.__define_notation()

    def __define_notation(self) -> ReactionNotation:
        if '(' in self.reaction and ')' in self.reaction:
            return ReactionNotation.SOVETIAN
        
        if '->' in self.reaction:
            return ReactionNotation.CHEMIST
        
        return ReactionNotation.UNDEFINED

    def to_reaction(self) -> Reaction:
        nucleus = self.split_nucleus()

        beam = self.to_nuclei(nucleus[0])
        target = self.to_nuclei(nucleus[1])
        fragment = self.to_nuclei(nucleus[2])

        return Reaction(beam, target, fragment, self.energy)

    def to_nuclei(self, name: str) -> Nuclei:
        nuclons = ReactionMaster.nuclons_from_name(name)
        charge = ReactionMaster.charge_from_name(name)
        return Nuclei(charge, nuclons)

    def split_nucleus(self) -> list[str]:
        no_spaces = self.reaction.replace(' ', '')

        if self.notation == ReactionNotation.UNDEFINED:
            raise ValueError('Reaction was written incorrectly.')

        if self.notation == ReactionNotation.SOVETIAN:
            return ReactionMaster.split_sovetian(no_spaces)

        if self.notation == ReactionNotation.CHEMIST:
            return ReactionMaster.split_chemistry(no_spaces)
    
    @staticmethod
    def split_sovetian(input: str) -> list[str]:
        halfs = input.split(',')

        first_half = halfs[0].split('(')
        second_half = halfs[1].split(')')

        beam, target = first_half[1], first_half[0]
        fragment, residual = second_half[0], second_half[1]

        return [beam, target, fragment, residual]
    
    @staticmethod
    def split_chemistry(input: str) -> list[str]:
        halfs = input.split('->')

        first_half = halfs[0].split('+')
        second_half = halfs[1].split('+')

        beam, target = first_half[1], first_half[0]
        fragment, residual = second_half[1], second_half[0]

        return [beam, target, fragment, residual]

    @staticmethod
    def nuclons_from_name(name: str) -> int:
        name = name.lower()
        if name in ['p', 'd', 't']:
            return ['p', 'd', 't'].index(name) + 1

        pretend = ''
        for i in name:
            if i.isdigit():
                pretend += i

        return int(pretend)

    @staticmethod
    def charge_from_name(name: str) -> int:
        name = ''.join([char.lower() for char in name if char.isalpha()])
        return NAME2CHARGE[name]


if __name__ == '__main__':
    master = ReactionMaster("Li7 + d -> Li6 + t", 14.5)

    print(master.to_reaction().residual)
