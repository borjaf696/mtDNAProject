# Borja :)
import os, subprocess
import altair as alt
import pandas as pd
import numpy as np
import shutil
from shutil import copyfile

class Utils:
    @staticmethod
    def cpfile(source, dest):
        copyfile(source, dest)
        return dest

    @staticmethod
    def getcwd():
        return os.getcwd()

    @staticmethod
    def chdir(dir):
        os.chdir(dir)

    @staticmethod
    def mkdir(dir):
        if not os.path.exists(dir):
            os.mkdir(dir)
        else:
            print('Directory: ',dir,' already exists')

    @staticmethod
    def exist_dir(d):
        return os.path.isdir(d)

    @staticmethod
    def exists(d):
        return os.path.exists(d)

    @staticmethod
    def remove_file(file):
        print('Removing: ',file)
        if Utils.exists(file):
            os.remove(file)

    @staticmethod
    def remove_dir(folder):
        shutil.rmtree(folder, ignore_errors=True)

    @staticmethod
    def get_dirs(path, bound = None):
        if bound is None:
            return [path+'/'+dI for dI in os.listdir(path) if os.path.isdir(os.path.join(path,dI))]
        else:
            dirs = []
            for dI in os.listdir(path):
                if len(dirs) == bound:
                    return dirs
                if os.path.isdir(os.path.join(path, dI)):
                    dirs.append(path+'/'+dI)
            return dirs

    @staticmethod
    def get_files(path, extension = ['fq','fastq','fasta','fa'], content = None):
        if content is None:
            return [path+'/'+t for t in os.listdir(path) if not os.path.isdir(path+'/'+t) and t.split('.')[-1] in extension]
        else:
            return [path+'/'+t for t in os.listdir(path) if not os.path.isdir(path+'/'+t) and content in t and t.split('.')[-1] in extension]

    @staticmethod
    def get_files_recursive(path, extension = ['fq','fa','fasta','fastq'], threshold = None):
        paths, results = [path], []
        while len(paths) > 0:
            p = paths[0]
            paths = paths[1:len(paths)]
            for t in os.listdir(p):
                n_path = p+'/'+t
                if os.path.isdir(n_path):
                    paths.append(n_path)
                elif t.split('.')[-1] in extension:
                        results.append(n_path)
                if threshold is not None and len(results) > threshold:
                    return results
        return results

    @staticmethod
    def get_files_recursive_content(path, content, threshold = None, avoid = None):
        paths, results = [path], []
        while len(paths) > 0:
            p = paths[0]
            paths = paths[1:len(paths)]
            for t in os.listdir(p):
                n_path = p + '/' + t
                if os.path.isdir(n_path):
                    paths.append(n_path)
                elif content in n_path and avoid is not None and avoid not in n_path:
                    results.append(n_path)
                elif content in n_path and avoid is None:
                    results.append(n_path)
            if threshold is not None and len(results) > threshold:
                return results
        return results

    @staticmethod
    def append_files(list_files, output_file):
        outDir = '/'.join(output_file.split('/')[:-1])
        if not Utils.exists(outDir):
            Utils.mkdir(outDir)
        with open(output_file,'w+') as out_file:
            for f in list_files:
                with open(f,'r') as f_read:
                    for line in f_read.readlines():
                        out_file.write(line)
        return output_file

    @staticmethod
    def append_files_bash(list_files, output_file):
        outDir = '/'.join(output_file.split('/')[:-1])
        if not Utils.exists(outDir):
            Utils.mkdir(outDir)
        cmd = ['cat']+[t for t in list_files]
        Utils.executecmd(cmd, output_file)
        return output_file

    @staticmethod
    def executecmd(args, out = None):
        '''
        :param path: LIST of args
        :return:
        '''
        print(' '.join(args))
        if out is None:
            subprocess.call(args)
        else:
            with open(out, 'w+') as fpout:
                subprocess.call(args, stdout = fpout)

    @staticmethod
    def write_list(lst,output, type = 'int'):
        with open(output, 'w+') as f:
            for l in lst:
                if type != 'int':
                    line = ','.join(map(lambda x:str(x),l))
                    f.write(line+'\n')
                else:
                    f.write(str(l)+'\n')

    @staticmethod
    def export_dict(dictionary,name, format = 'csv'):
        if format == 'csv':
            with open(name, 'w+') as f:
                for key, val in dictionary.items():
                    f.write(str(key))
                    for k,v in val.items():
                        f.write(','+str(v))
                    f.write('\n')
        if format == 'fasta':
            with open(name, 'w+') as f:
                for key, val in dictionary.items():
                    f.write(str(key))
                    f.write(str(val))

    @staticmethod
    def writeMatrix(m, file = None):
        assert file is not None

        rows,cols = np.shape(m)
        with open(file, 'w+') as f:
            for i in range(rows):
                for j in range(cols):
                    f.write(str(m[i,j])+' ')
                f.write('\n')


class StatsReport:
    @staticmethod
    def show_histogram(histogram, path, y_label = 'y', x_label = 'x'):
        d = {y_label:histogram,x_label:range(len(histogram))}
        source = pd.DataFrame(d)
        chart = alt.Chart(source).mark_bar().encode(alt.X(x_label, bin = True),y=y_label)
        chart_2 = alt.Chart(source).mark_bar().encode(alt.X(x_label,  bin=alt.Bin(maxbins=100)), y = y_label)
        chart.save(path)
        chart_2.save(path+'_100_bins.html')

    @staticmethod
    def dotPlot(df, X, Y, nameFile, color = None):
        chart = alt.Chart(df)
        chart = chart.mark_point().encode(x=X, y = Y) if color is None else chart.mark_point().encode(x=X,y=Y, color=color)
        chart.save(nameFile+'.html')

    @staticmethod
    def heatmap(df, X, Y, m, nameFile):
        chart = alt.Chart(df).mark_rect().encode(x=X+':O',y=Y+':O',color=m+':Q')
        chart.save(nameFile+'.html')

    @staticmethod
    def iterative_levenshtein(s, t, costs=(1, 1, 1)):
        rows = len(s) + 1
        cols = len(t) + 1
        deletes, inserts, substitutes = costs

        dist = [[0 for x in range(cols)] for x in range(rows)]
        for row in range(1, rows):
            dist[row][0] = row * deletes
        for col in range(1, cols):
            dist[0][col] = col * inserts

        for col in range(1, cols):
            for row in range(1, rows):
                if s[row - 1] == t[col - 1]:
                    cost = 0
                else:
                    cost = substitutes
                dist[row][col] = min(dist[row - 1][col] + deletes,
                                     dist[row][col - 1] + inserts,
                                     dist[row - 1][col - 1] + cost)
        '''for r in range(rows):
            print(dist[r])'''
        return dist[rows-1][cols-1]
