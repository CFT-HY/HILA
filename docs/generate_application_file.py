import os,sys,re

def get_source_files(file_path):
    source_files = []
    for root, dirs, files in os.walk(file_path):
        if 'src' in root:
            for file in files:
                source_files.append(root + '/' + file)
                
    return source_files

def find_brief(source_files):
    files_and_briefs = {}
    for file in source_files:
        opened_file = open(file,"r")
        file_content = opened_file.read()
        if re.search(r'@file', file_content):
            print('file is documented')
            start = int(re.search(r'@brief',file_content).span()[1])
            end = int(re.search(r'@details', file_content).span()[0])
            brief = repr(file_content[start:end])
            brief = re.sub(r'\\n \*', '',brief)
            #remove first and last character ' with indexing
            files_and_briefs[file] = str(brief[1:-1])
        
    return files_and_briefs
    
def generate_application_file(files_and_briefs,write_file):
    f = open(write_file,'w')
    f.write("/** \n\n \page hila_applications HILA Applications \n\n__List of documented applications/simulations offered by HILA.__\n\n")
    for file,brief in files_and_briefs.items():
        f.write("{}:\n\n".format(file.split('/')[-1]))
        f.write("> {}\n\n".format(brief))
    f.write("*/")
    f.close()

            
        
if __name__ == "__main__":
    root = sys.argv[1]
    write_file = sys.argv[2]
    source_files = get_source_files(root)
    files_and_briefs = find_brief(source_files)
    generate_application_file(files_and_briefs,write_file)
    
    