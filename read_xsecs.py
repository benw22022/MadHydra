import glob
import os

def main():
    
    logfiles = glob.glob("multirun/*/*/log.generate")
    
    for fpath in logfiles:
        with open(fpath, 'r') as file:
            for line in file:
                if "Cross-section :" in line:
                    print(os.path.basename(os.path.dirname(fpath)), line)
                    break
            print(os.path.basename(os.path.dirname(fpath)), "     NO CROSS-SECTION FOUND")
    
if __name__ == "__main__":
    
    main()