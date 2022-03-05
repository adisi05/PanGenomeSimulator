from multiprocessing import Process, Manager, Pool


def test(num, q):
    q.put(num)


if __name__ == '__main__':
    # https://stackoverflow.com/questions/51398008/pool-and-manager-in-multiprocessing-module
    pool = Pool(4)
    q = Manager().queue()
    for i in range(30):
        pool.apply_async(test, (i, q))
    pool.close()
    pool.join()