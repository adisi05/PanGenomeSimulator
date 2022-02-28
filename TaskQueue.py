from multiprocessing import Process
from multiprocessing import Queue
import time
# from messages import MessageStream
# from utils.decode_and_encode import decode_message
#
# FULL_NAME: str = 'Adi Sivan'
# PHONE_NUMBER: str = '0526053358'
# PERSONAL_COMMENT: str = 'I know I have a problem deleting messages in line 89, the published_counter is not correct when the last sequence was published and deleted. I didnt have time for that unfortunately :(' \
#                         'Besides that the exercise was OK. I think more time could be useful'
#
# def consume_process_and_produce(stream: MessageStream, sequence_min_length: int = 0):
#     """
#     Use stream: MessageStream argument to consume and produce processed messages.
#
#     Step1 (decode + speed + order):
#         The messages are base64 encoded strings and MUST be:
#          - decoded before published
#          - produced in the same order they were consumed
#          - published as early as possible
#
#         for example, for listen input of: 1,2,3
#         output should be: 1,2,3
#
#     Step2 (Step1 + sequence lockout):
#         - Keep all constraints from Step1
#         - Publish only sequences with len >= 3 (example). Sequence is consecutive repetition of the same message
#
#         for example, for listen input of: 1,4,2,2,5,5,5,9,9,8,8,8,8,6
#         output should be: 5,5,5,8,8,8,8
#
#     """
#
#     threads_count = 5
#     task_queue = TaskQueue(decode_message, threads_count)
#     task_counter = 0
#     published_counter = 0
#
#     # Step 1
#     if sequence_min_length == 0:
#         for m in stream.listen():
#             if m:
#                 task_counter += 1
#                 task_queue.add_task(m,task_counter)
#             published_counter = publish_single_message(published_counter, stream, task_queue)
#         publish_remaining_single_messages(published_counter, stream, task_counter, task_queue)
#
#     # Step 2
#     else:
#         sequence_counter = 0
#         last_message = -1
#         for m in stream.listen():
#             if m:
#                 task_counter += 1
#                 task_queue.add_task(m, task_counter)
#             published_counter, sequence_counter, last_message = publish_sequence_message(published_counter, stream, task_queue, sequence_counter, last_message, sequence_min_length)
#         publish_remaining_sequence_messages(published_counter, stream, task_queue, sequence_counter, last_message, sequence_min_length, task_counter)
#
# def publish_remaining_single_messages(published_counter, stream, task_counter, task_queue):
#     while published_counter < task_counter:
#         published_counter = publish_single_message(published_counter, stream, task_queue)
#
# def publish_single_message(published_counter, stream, task_queue):
#     if published_counter + 1 in task_queue._decoded_messages:
#         published_counter += 1
#         stream.publish(task_queue._decoded_messages[published_counter])
#         del task_queue._decoded_messages[published_counter]
#     return published_counter
#
# def publish_remaining_sequence_messages(published_counter, stream, task_queue, sequence_counter, last_message, sequence_min_length, task_counter):
#     while published_counter < task_counter:
#         published_counter, sequence_counter, last_message = publish_sequence_message(published_counter, stream, task_queue, sequence_counter, last_message, sequence_min_length)
#
# def publish_sequence_message(published_counter, stream, task_queue, sequence_counter, last_message, sequence_min_length):
#     if published_counter + 1 in task_queue._decoded_messages:
#         m = task_queue._decoded_messages[published_counter+1]
#         if m.decoded_data == last_message:
#             sequence_counter += 1
#             if sequence_counter == sequence_min_length:
#                 # publish last messages
#                 for i in range(sequence_min_length):
#                     published_counter = publish_single_message(published_counter, stream, task_queue)
#             elif sequence_counter > sequence_min_length:
#                 # publish only last
#                 published_counter = publish_single_message(published_counter, stream, task_queue)
#         else:
#             last_message = m.decoded_data
#             # TODO delete unpublished messages, currently only ignoring them...
#             published_counter += sequence_counter
#             sequence_counter = 1
#     return published_counter, sequence_counter, last_message
#
# if __name__ == "__main__":
#     from test.test_runner import run_all_tests, registration
#
#     registration(FULL_NAME, PHONE_NUMBER, PERSONAL_COMMENT)
#
#     run_all_tests()

# Reference: https://medium.com/@shashwat_ds/a-tiny-multi-threaded-job-queue-in-30-lines-of-python-a344c3f3f7f0
class TaskQueue(Queue):

    def __init__(self, func, num_workers=1):
        Queue.__init__(self)
        self._num_workers = num_workers
        self.start_workers()
        self._func = func

    def add_task(self, *arg):
        self.put(*arg)

    def start_workers(self):
        for i in range(self._num_workers):
            t = Process(target=self.worker)
            t.daemon = True
            t.start()

    def worker(self):
        while True:
            self._func(self.get())
            self.task_done()
