import multiprocessing


class Logger:

    def __init__(self, is_concurrent: bool = False, debug: bool = False):
        self.is_concurrent = is_concurrent
        self.debug = debug

    def set_concurrent(self, is_concurrent):
        self.is_concurrent = is_concurrent

    def message(self, message: str):
        log_message(message, self.is_concurrent)

    def debug_message(self, message: str):
        if self.debug:
            log_message(message, self.is_concurrent)


def log_message(message: str, is_concurrent: bool = False):
    if is_concurrent:
        curr_proc = multiprocessing.current_process()
        print(curr_proc.name, ' : ', message)
    else:
        print(message)
