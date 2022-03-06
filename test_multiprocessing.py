########################################## Old stuff: ##########################################

# import time
# from multiprocessing import Process, Manager, Pool, Queue
# from multiprocessing.managers import SyncManager
# from queue import PriorityQueue
#
# # https://stackoverflow.com/questions/53613583/queue-function-for-multiprocess-priority-queue-in-python-with-syncmanager-class
# # Python 3.7.1
# from queue import PriorityQueue
# from multiprocessing.managers import SyncManager
#
#
# SENTINEL = None
#
#
# class MyPriorityQueue(PriorityQueue):
#     def get_attribute(self, name):
#         return getattr(self, name)
#
#
# class MyManager(SyncManager):
#     pass
#
#
# def get_manager():
#     MyManager.register("PriorityQueue", MyPriorityQueue)
#     m = MyManager()
#     m.start()
#     return m
#
# # def f(q):
# #     for item in iter(lambda: q.get()[1], SENTINEL):
# #         print(item)
# #     print(f'queue: {q.get_attribute("queue")}')
#
#
#
# def f(q):
#     # while not q.get()[1].empty():
#     #     item = q.get()[1].get()
#     #     if item is None:
#     #         break
#     #     print("{} removed {} from the queue".format(pool.current_thread(), item))
#     #     # queue.task_done()
#     #     time.sleep(1)
#
#     for item in iter(lambda: q.get()[1], SENTINEL):
#         print(item)
#     print(f'queue: {q.get_attribute("queue")}')
#
# def f2(*args):
#     print(args)
#     print(args[0]+args[1]+args[2])
#
#
# if __name__ == '__main__':
#
#     m = get_manager()
#     pq = m.PriorityQueue()
#
#     tasks = enumerate([f'item_{i}' for i in range(5)] + [SENTINEL])
#
#     for task in tasks:
#         pq.put(task)
#
#     print(f'queue: {pq.get_attribute("queue")}')
#     print(f'maxsize: {pq.get_attribute("maxsize")}')
#
#     # https://stackoverflow.com/questions/51398008/pool-and-manager-in-multiprocessing-module
#     # pool = Pool(4)
#     # for i in range(30):
#     #     pool.apply_async(f,args=(pq,))
#
#     # https://github.com/MayroseLab/chimeraBuster/blob/master/detect_chimeric_genes.py
#     ncpu = 4
#     with Pool(ncpu) as pool:
#         r = pool.apply_async(f2,args=("test","hello","world"))
#
#     # pool.close()
#     # pool.join()


########################################## First example: ##########################################
# https://stackoverflow.com/questions/35908987/multiprocessing-map-vs-map-async

from multiprocessing import Pool

def f(x):
    time.sleep(1)
    print(x*x) #TODO hsould acquire lock

def main1():
    pool = Pool(processes=4)
    # pool.map(f, range(10))
    r = pool.map_async(f, range(10))
    # DO STUFF
    print('HERE')
    print('MORE')
    r.wait()
    print('DONE')

#--------------------------------------------------------------------------------#

########################################## Second example: ##########################################
# https://stackoverflow.com/a/3232026
import multiprocessing, time

def task(args):
    count = args[0]
    queue = args[1]
    for i in range(count):
        queue.put("%d mississippi" % i)
    return "Done"


def main2():
    manager = multiprocessing.Manager()
    q = manager.Queue()
    pool = multiprocessing.Pool()
    result = pool.map_async(task, [(x, q) for x in range(10)])
    time.sleep(1)
    while not q.empty():
        print(q.get())
    print(result.get())

#--------------------------------------------------------------------------------#

if __name__ == "__main__":
    main2()
