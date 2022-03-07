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

def my_task1(x):
    time.sleep(1)
    print(x*x) #TODO hsould acquire lock

def main1():
    pool = Pool(processes=4)
    # pool.map(f, range(10))
    r = pool.map_async(my_task1, range(10))
    # DO STUFF
    print('HERE')
    print('MORE')
    r.wait()
    print('DONE')

########################################## Second example: ##########################################
# https://stackoverflow.com/a/3232026

import multiprocessing, time

def my_task2(args):
    count = args[0]
    queue = args[1]
    for i in range(count):
        queue.put("%d mississippi" % i)
    return "Done"


def main2():
    manager = multiprocessing.Manager()
    q = manager.Queue()
    pool = multiprocessing.Pool()
    result = pool.map_async(my_task2, [(x, q) for x in range(10)])
    time.sleep(1)
    while not q.empty():
        print(q.get())
    print(result.get())

########################################## Third example: ##########################################

def my_task3(args):
    x = args[0]
    queue = args[1]
    queue.put("%d mississippi" % x)
    while not queue.empty():
        print(queue.get())
        time.sleep(5)
    return "Done"

def main3():
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(processes=4)

    q = manager.Queue()
    half_amount = 5
    for i in range(half_amount):
        q.put("task number %d" % i)
    result = pool.map_async(my_task3, [(x, q) for x in range(10)])
    time.sleep(5)
    for i in range(half_amount, 2*half_amount):
        q.put("task number %d" % i)
    result.wait()
    print(result.get())

########################################## Forth example: ##########################################

def my_task4(x, queue):
    curr_proc = multiprocessing.current_process()
    print('Hi there! current process:', curr_proc.name, curr_proc._identity)
    queue.put("%d mississippi" % x)
    while not queue.empty():
        print(queue.get())
        time.sleep(5)
    return "Done"

def main4():
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(processes=4)

    q = manager.Queue()
    half_amount = 5
    re = [0]*10
    for i in range(half_amount):
        q.put("task number %d" % i)
    for x in range(10):
        re[x]=pool.apply_async(my_task4, args=(x, q))
    for i in range(half_amount, 2*half_amount):
        q.put("task number %d" % i)
    pool.close()
    pool.join()



########################################## Fifth example: ##########################################

def my_task5(lock, x, queue):


    curr_proc = multiprocessing.current_process()
    print('Hi there! current process:', curr_proc.name, curr_proc._identity)
    queue.put("%d mississippi" % x)
    bye = False
    while not queue.empty() and not bye:
        lock.acquire()
        if queue.empty():
            bye = True
        else:
            print(queue.get())
        lock.release()
        if bye:
            break
        time.sleep(5)
    return "Done"

def main5():
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(processes=4)
    lock = manager.Lock()
    q = manager.Queue()
    half_amount = 5
    re = [0]*10
    for i in range(half_amount):
        q.put("task number %d" % i)
    for x in range(10):
        pool.apply_async(my_task5, args=(lock, x, q))
    for i in range(half_amount, 2*half_amount):
        q.put("task number %d" % i)
    pool.close()
    pool.join()

# def main6():
#     manager = multiprocessing.Manager()
#     # pool = multiprocessing.Pool(processes=4)
#     lock = manager.Lock()
#     q = manager.Queue()
#     half_amount = 5
#     ncpu = 4
#     with multiprocessing.Pool(processes=ncpu) as pool:
#         for i in range(half_amount):
#             q.put("task number %d" % i)
#         for x in range(10):
#             pool.apply_async(my_task5, args=(lock, x, q))
#         for i in range(half_amount, 2 * half_amount):
#             q.put("task number %d" % i)
#     pool.close()
#     pool.join()


if __name__ == "__main__":
    main5()
