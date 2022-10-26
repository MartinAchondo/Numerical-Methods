import subprocess
import time


start_time = time.time()
script_name = "main"
subprocess.run(['matlab','-batch',script_name])
print("--- %s seconds ---" % (time.time() - start_time))