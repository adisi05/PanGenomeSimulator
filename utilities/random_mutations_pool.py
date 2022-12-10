import random
from utilities.common_data_structues import MutType, Region
from utilities.logger import Logger


class RandomMutationsPool:
    def __init__(self, indels_per_region: dict, snps_per_region: dict, max_mutations_in_window: int,
                 sv_list: list = [], logger: Logger = None):
        self.logger = logger if logger else Logger()
        # TODO add SVs
        self.indels_per_region = indels_per_region
        self.snps_per_region = snps_per_region
        self.options = {}
        for region_name, count in indels_per_region.items():
            if count != 0:
                self.options[(MutType.INDEL.value, region_name)] = count
        for region_name, count in snps_per_region.items():
            if count != 0:
                self.options[(MutType.SNP.value, region_name)] = count
        self.overall_count = min(max_mutations_in_window, sum(self.options.values()))

        self.logger.debug_message(f"Created random mutations pool")
        self.logger.debug_message(f"Overall planned random mutations within window: {self.overall_count}")
        self.logger.debug_message(f"Mutations distribution: {list(self.options.items())}")

    def has_next(self) -> bool:
        return self.overall_count > 0

    def get_next(self) -> (MutType, Region):
        if self.overall_count <= 0:
            return None
        option_with_count_list = self.options.items()
        # https://pynative.com/python-weighted-random-choices-with-probability/
        choice = random.choices([opt[0] for opt in option_with_count_list],
                                weights=[opt[1] for opt in option_with_count_list],
                                k=1)[0]
        self.overall_count -= 1
        self.options[choice] -= 1
        if self.options[choice] == 0:
            del self.options[choice]
        mut_type_name, region_name = choice

        self.logger.debug_message(f"Remained mutations count within window: {self.overall_count}")

        return MutType(mut_type_name), Region(region_name)
