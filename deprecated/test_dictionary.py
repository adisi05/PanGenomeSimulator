import multiprocessing
from queue import Queue
from time import sleep

def func(dict):
    dict_a = {}
    dict_b = {}
    for k in sorted(dict.keys()):
        dict[k].sort()
        l = len(dict[k])
        mid = round(l/2)
        dict_a[k] = dict[k][0:min(mid,l)]
        dict_b[k] = dict[k][min(mid,l):l]
    return dict, dict_a, dict_b

def print_dict(d):
    sleep(10)
    print("Dictionary:",d)

if __name__ == "__main__":
    dict1 = {'b': ['b1', 'b3', 'b2'],
            'c': ['c3','c1','c2'],
            'a': ['a1','a3'],
            'd': ['d4','d2']}

    print("1", dict1)
    dict2, dict2_a, dict2_b = func(dict1)
    dict2['b'].append('b4')
    print("2", dict1)
    print("3", dict2)
    print("4", dict2_a)
    print("5", dict2_b)
    q = Queue()
    q.put(dict1)
    q.put(dict2_a)
    q.put(dict2_b)
    pool = multiprocessing.Pool(processes=3)
    for d in q.queue:
        pool.apply_async(print_dict, args=(d,))
    dict2_a['e'] = ['e5', 'e6']
    pool.close()
    pool.join()
