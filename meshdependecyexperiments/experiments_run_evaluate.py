import meshio
import sys
from pylab import *
import matplotlib as mpl
import os
import time
import numpy

def open_out_dir(out_dir):
    x_max = []
    y_max = []
    z_max = []
    x_min = []
    y_min = []
    z_min = []
    
    for i in range(1):    
        mesh_i = meshio.read(out_dir+'file_'+str(i)+'.vtu')
        #converte para array de numpy e calcula a transversa, para que a coluna represente os eixos x, y e z.
        arr_pnts = numpy.transpose(numpy.array(mesh_i.points))
        #calcular os mínimos de cada eixo
        #calcular o máximo de cada eixo
def get_output_dir(param_file):
    f = open(param_file, "r")
    if f:
        for linha in f:
            pass
        f.close()
        return linha;
    else:
        print("Param file can not be opened: "+param_file)
def run_experiment(param_file):
    exe_path   = 'D:\\Ricardo\\OneDrive\\Documentos\\GitHub\\FisioPacerGPU\\FisioPacerGPU\\x64\Release\\FisioPacerGPU.exe'
    start = time.time()
    stream = os.popen(exe_path+' '+param_file)
    output = stream.read()
    print(output)
    end = time.time()
    run_time_s = end - start
    print('Tempo: '+ str(run_time_s)+'s')
    out_dir = get_output_dir(param_file)
    open_out_dir(out_dir)
def run_experiments():
    
    # Using system() method to
    # execute shell commands
    
    param_exp1 = 'D:\\Ricardo\\OneDrive\\Documentos\\GitHub\\FisioPacerGPU\\meshdependecyexperiments\\cubo74\\exp_gravidade.param'
    run_experiment(param_exp1)

    param_exp2 = 'D:\\Ricardo\\OneDrive\\Documentos\\GitHub\\FisioPacerGPU\\meshdependecyexperiments\\cubo295\\exp_gravidade.param'
    #run_experiment(param_exp2)

    
run_experiments();

#meshG = meshio.read('grossa.vtu', 'vtu-ascii')
#meshR = meshio.read('refinada.vtu', 'vtu-ascii')
#print('Ponto:')
#data = {'x':[len(meshG.points), len(meshR.points)],
#            'y':[pontoExp1(meshG.points), pontoExp1(meshR.points)]}
#print(data)
