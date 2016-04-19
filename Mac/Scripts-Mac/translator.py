import subprocess, platform, os,shlex




if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")



def translate(fileName):
    shell = "perl ~/Desktop/bioinformatics/Scripts/translatoX.pl "
    shell = shell + "-i " + fileName
    shell = shell + " -o " + fileName
    shell = shell + ' -t T -c 1 -g "-b5 a"'
    print(shell)
    os.system(shell)

    return


def crawl(folder):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if "_rn.fasta" in fileName and ".DS_Store" not in fileName and "clean" not in fileName:
                translate(fileName)

    return


crawl("./project_temp/dsDNA/Apapillomavirus/NP_041332.1")
