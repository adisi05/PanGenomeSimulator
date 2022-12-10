import multiprocessing


class Logger:

    def __init__(self, is_concurrent=False):
        self.is_concurrent = is_concurrent

    def set_concurrent(self, is_concurrent):
        self.is_concurrent = is_concurrent

    def log(self, message: str):
        if self.is_concurrent:
            curr_proc = multiprocessing.current_process()
            print(curr_proc.name, ' : ', message)
        else:
            print(message)
