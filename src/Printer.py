from timeit import default_timer as timer
import time

class Printer():
    def __init__(self):
        self.startTime = timer()

    def __call__(self, printText):
        readableCurrentTime = self.determine_time_since_start_in_human_readable_format(timer())
        print("\n----> " + str(readableCurrentTime) + ": " + printText + "\n", flush=True)

    def determine_time_since_start_in_human_readable_format(self, currentTime):
        return time.strftime("%H:%M:%S", time.gmtime(currentTime- self.startTime))


