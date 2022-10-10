class SimulationWindow:
    def __init__(self, debug: bool = False):

        self.windows_start_offset: int = 0
        self.start: int = 0
        self.end: int = 0
        self.original_end: int = 0
        self.debug = debug

    def next_window(self, new_start: int = -1, new_end: int = -1):
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
        if self.debug:
            print(f"Updated window. Overall windows offset is {self.windows_start_offset}"
                  f" and therefore the adjusted parameters are window: start={self.start}, end={self.end}")

    def adjust_window(self, end_shift: int = 0):
        if self.debug:
            print(f"Adjusting window, end shift is {end_shift}")
        self.end += end_shift
