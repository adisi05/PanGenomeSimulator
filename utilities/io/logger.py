import multiprocessing


class Logger:

    def __init__(self, name: str = '', debug: bool = False):
        self.name = name
        self.debug = debug

    def message(self, message: str):
        log_message(message, self.name)

    def debug_message(self, message: str):
        if self.debug:
            log_message(message, self.name)


def log_message(message: str, name: str = ''):
    if len(name) > 0:
        print(name, ' : ', message)
    else:
        print(message)
