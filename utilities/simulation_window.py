from utilities.common_data_structues import Mutation, MutType


class SimulationWindow:
    def __init__(self, debug: bool = False):

        self.windows_start_offset: int = 0
        self.start: int = 0
        self.end: int = 0
        self.original_end: int = 0
        self.debug = debug
        # Blocklist explanation:
        # blocklist[pos] = 0		safe to insert variant here
        # blocklist[pos] = 1		indel inserted here
        # blocklist[pos] = 2		snp inserted here
        # blocklist[pos] = 3		invalid position for various processing reasons
        self.blocklist = {}  # TODO re-consider !!!

    def next_window(self, new_start: int = -1, new_end: int = -1):
        # TODO sanity check about start, end. Maybe consider N regions?
        last_window_offset = self.end - self.original_end
        new_start += (self.windows_start_offset + last_window_offset)
        new_end += (self.windows_start_offset + last_window_offset)
        if new_start < 0 or new_end < 0 or not (self.end <= new_start < new_end):
            raise Exception(f'Illegal window positions. Start and end should be positive integers. '
                            f'The new window should come after the previous one and not overlap with it. '
                            f'Got the next values: start={new_start}, end={new_end}, '
                            f'while the previous values are: start={self.start}, end={self.end}')
        self.windows_start_offset += last_window_offset
        self.start = new_start
        self.end = new_end
        self.original_end = self.end
        self.blocklist = {}  # np.zeros(end-start, dtype='<i4')
        if self.debug:
            print(f"Updated window. Overall windows offset is {self.windows_start_offset}"
                  f" and therefore the adjusted parameters are window: start={self.start}, end={self.end}")

    def adjust_window(self, end_shift: int = 0):
        # TODO sanity check about start, end. Maybe consider N regions?
        if self.debug:
            print(f"Adjusting window, end shift is {end_shift}")
        self.end += end_shift

    def update_blocklist(self, mutation: Mutation):
        for k in range(mutation.position, mutation.position + len(mutation.ref_nucl)):  # TODO continue here
            self.blocklist[k] = 1 if mutation.mut_type.value == MutType.INDEL.value else 2

    def check_blocklist(self, mutation: Mutation):
        for k in range(mutation.position, mutation.position + len(mutation.ref_nucl)):
            if self.blocklist.get(k, 0):
                return False
        return True