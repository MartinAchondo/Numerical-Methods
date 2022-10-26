import time
import matlab.engine

def take_time(f):
    def wrapper(*kargs,**kwargs):
        start_time = time.time()
        f()
        print("--- %s seconds ---" % (time.time() - start_time))
    return wrapper

@take_time
def main():
    eng = matlab.engine.start_matlab()
    eng.main(nargout=0)
    eng.quit()

if __name__=='__main__':
    main()