
# Você pode encontrar informações sobre este trabalho no seguinte repositório: https://github.com/luandebrito/PythonGUI_for_ProcessControl
# Se possível, cite os trabalhos encontrados la no seu projeto :)

# You can find information about this work in the following repository: https://github.com/luandebrito/PythonGUI_for_ProcessControl
# If possible, cite the works found there in your project :)


#importar bibliotecas e pacotes // import libraries and packages
from PyQt5 import QtWidgets, QtCore, uic
from PyQt5 import uic,QtGui
from PyQt5.QtCore import *
import pyqtgraph as pg
import sys
import nidaqmx
from nidaqmx.constants import TerminalConfiguration
import time
from xlwt import Workbook
from datetime import datetime
import numpy as np
from scipy import integrate, interpolate, optimize, signal
from scipy.ndimage import gaussian_filter1d

#Thread de sintonia do controlador // Controller tuning thread
class TuningThread(QThread): 

    change_value = pyqtSignal(list)
    change_value2 = pyqtSignal(int)
    chave = True
    var = 0
    pv = []
    output = []
    tdegrau = 0
    y0 = 0
    count = 0
    tlist = []
    t1 = 0
    t2 = 0
    Adegrau = 0
    kp = 0
    taup = 0
    thetap = 0
    A1 = 0
    Asteady = 0
    A0 = 0
    A2 = 0
    d2 = []
    getD1 = False
    alfa = 0.2

    def run(self):
         
        TuningThread.pv.clear()
        TuningThread.output.clear()
        TuningThread.tlist.clear()
        TuningThread.count = 0
    
        t = 0    

        while TuningThread.chave == True:
          
            TuningThread.pv.append(WorkerThread1.nivel[WorkerThread1.indice])
            TuningThread.output.append(WorkerThread1.bomba[WorkerThread1.indice])
            TuningThread.tlist.append(t) 
            self.change_value.emit(TuningThread.pv)
            t = t+ WorkerThread1.dt
            TuningThread.count = TuningThread.count + 1
            time.sleep(WorkerThread1.dt)


    def getData1(self):
        TuningThread.getD1 = True 
        
    def getData2(self):
            
        TuningThread.t2 = TuningThread.count - 1
        data = TuningThread.pv[TuningThread.t1:TuningThread.t2+1]
        t = TuningThread.tlist[TuningThread.t1:TuningThread.t2+1]
        u = TuningThread.output[TuningThread.t1:TuningThread.t2+1]
        dataexp = [x - TuningThread.pv[TuningThread.t1] for x in data]
        TuningThread.kp = (TuningThread.pv[TuningThread.t2]-TuningThread.pv[TuningThread.t1]) / TuningThread.Adegrau
            
        if (len(dataexp) % 2) == 0:
            l = len(dataexp)-1
        else:
            l = (len(dataexp))

        smooth_data = signal.savgol_filter(dataexp, 15, 1)
        d1 = np.gradient(smooth_data, WorkerThread1.dt) #primeira derivada // first derivative
        TuningThread.d2 = np.gradient(d1, WorkerThread1.dt) #segunda derivada // second derivative
        
        try:
            infls = np.where(np.gradient(np.sign(TuningThread.d2),WorkerThread1.dt))[0]
        except:
            infls = [100]
                
        f = 0
        f3 = interpolate.splrep(t,dataexp, s = 15, k = 3)  # se necesario (if its necessary) fill_value = "extrapolate"
        fy = interpolate.splev(t , f3)

        d_while = True
    
        while d_while == True:
            if d1[infls[f]]>0 and TuningThread.Adegrau>0:
                d_while = False
            elif d1[infls[f]]<0 and TuningThread.Adegrau<0:
                d_while = False
            else:
                f = f + 1
    
        TuningThread.thetap = (t[infls[f+1]] - (dataexp[infls[f+1]] / d1[infls[f+1]])) - t[0]
        if TuningThread.thetap < 0:
            TuningThread.thetap = 0

        z = 0.632 * (TuningThread.pv[TuningThread.t2]-TuningThread.pv[TuningThread.t1])

        absolute_difference_function = lambda list_value : abs(list_value - z)
        zvalue = min(dataexp, key=absolute_difference_function)
        zindex = dataexp.index(zvalue)
        TuningThread.taup = t[zindex] - TuningThread.thetap - t[0]

        u = TuningThread.output
        u0 = u[0]
        yp0 = TuningThread.pv[0]

        # specify number of steps // especificar número de passos
        ns = len(TuningThread.tlist)
        delta_t = t[1]-t[0]
        # create linear interpolation of the u data versus time // criar interpolação linear
        uf = interpolate.interp1d(TuningThread.tlist,u)

        # define first-order plus dead-time approximation   // definir aproximação do modelo FOPDT  
        def fopdt(y,t,uf,Km,taum,thetam):
            # arguments
            #  y      = output // saída
            #  t      = time // tempo
            #  uf     = input linear function (for time shift) // função linear de entrada
            #  Km     = model gain // ganho do processo ou modelo
            #  taum   = model time constant // constante de tempo do processo
            #  thetam = model time-delay constant // constante de atraso do processo
            # time-shift u
            try:
                if (t-thetam) <= 0:
                    um = uf(0.0)
                else:
                    um = uf(t-thetam)
            except:
                #print('Error with time extrapolation: ' + str(t))
                um = u0
            # calculate derivative // calcular derivada
            dydt = (-(y-yp0) + Km * (um-u0))/taum
            return dydt

        # simulate FOPDT model with x=[Km,taum,thetam] // simular modelo FOPDT
        def sim_model(x):
            # input arguments // argumentos de entrada
            Km = x[0]
            taum = x[1]
            thetam = x[2]
            # storage for model values // armazenamento de valores do modelo
            ym = np.zeros(ns)  # model // modelo
            # initial condition // condição inicial
            ym[0] = yp0
            # loop through time steps // loop entre passos de tempo
            for i in range(0,ns-1):
                ts = [TuningThread.tlist[i],TuningThread.tlist[i+1]]
                y1 = integrate.odeint(fopdt,ym[i],ts,args=(uf,Km,taum,thetam))
                ym[i+1] = y1[-1]
            return ym

        x1 = np.zeros(3)
        x1[0] = TuningThread.kp #2.0 # Km
        x1[1] = TuningThread.taup #3.0 # taum
        x1[2] = TuningThread.thetap # thetam
        ym1 = sim_model(x1)


        TuningWindow.dados.setData(TuningThread.tlist , TuningThread.pv)
        TuningWindow.guess.setData(TuningThread.tlist , ym1)
        TuningWindow.fopdt_optimized.clear()

        #if TuningThread.thetap <= 0:
        #    WorkerThread1.kc = 1/TuningThread.kp        #regras IMC para tempo morto = 0
        #    WorkerThread1.tauI = TuningThread.taup
        #    WorkerThread1.tauD = 0

        #else:  #imc agressivo 0,8Taup  moderado 8*taup  conservativo 80taup
        if TuningThread.thetap <= 0:
            WorkerThread1.kc = (1/TuningThread.kp)*((TuningThread.taup+0.5*TuningThread.thetap)/(0.1*TuningThread.taup+0.5*TuningThread.thetap))
            WorkerThread1.tauI = TuningThread.taup + 0.5*TuningThread.thetap
            WorkerThread1.tauD = (TuningThread.taup*TuningThread.thetap)/(2*TuningThread.taup+TuningThread.thetap)

        else:       #regras de ziegler-nichols
            WorkerThread1.kc = (1.2*TuningThread.taup)/(TuningThread.thetap*TuningThread.kp)
            WorkerThread1.tauI = 2*TuningThread.thetap
            WorkerThread1.tauD = TuningThread.thetap/2

        
        TuningWindow.gainp.setValue(TuningThread.kp)
        TuningWindow.taup.setValue(TuningThread.taup)
        TuningWindow.tetap.setValue(TuningThread.thetap)
        TuningWindow.gainc.setValue(WorkerThread1.kc)
        TuningWindow.tauc.setValue(WorkerThread1.tauI)
        TuningWindow.tetac.setValue(WorkerThread1.tauD)
            
    
    def update_model(self):
        TuningThread.t2 = TuningThread.count - 1
        data = TuningThread.pv[TuningThread.t1:TuningThread.t2+1]
        t = TuningThread.tlist[TuningThread.t1:TuningThread.t2+1]
        u = TuningThread.output[TuningThread.t1:TuningThread.t2+1]
        dataexp = [x - TuningThread.pv[TuningThread.t1] for x in data]
                    

        
        if TuningThread.thetap < 0:
            TuningThread.thetap = 0

        z = 0.632 * (TuningThread.pv[TuningThread.t2]-TuningThread.pv[TuningThread.t1])

        u = TuningThread.output
        u0 = u[0]
        yp = TuningThread.pv
        yp0 = TuningThread.pv[0]

        # specify number of steps
        ns = len(TuningThread.tlist)
        delta_t = t[1]-t[0]
        
        # create linear interpolation of the u data versus time
        uf = interpolate.interp1d(TuningThread.tlist,u)

        # define first-order plus dead-time approximation    
        def fopdt(y,t,uf,Km,taum,thetam):
            # arguments
            #  y      = output
            #  t      = time
            #  uf     = input linear function (for time shift)
            #  Km     = model gain
            #  taum   = model time constant
            #  thetam = model time constant
            # time-shift u
            try:
                if (t-thetam) <= 0:
                    um = uf(0.0)
                else:
                    um = uf(t-thetam)
            except:
                #print('Error with time extrapolation: ' + str(t))
                um = u0
            # calculate derivative
            dydt = (-(y-yp0) + Km * (um-u0))/taum
            return dydt

        # simulate FOPDT model with x=[Km,taum,thetam]
        def sim_model(x):
            # input arguments
            Km = x[0]
            taum = x[1]
            thetam = x[2]
            # storage for model values
            ym = np.zeros(ns)  # model
            # initial condition
            ym[0] = yp0
            # loop through time steps    
            for i in range(0,ns-1):
                ts = [TuningThread.tlist[i],TuningThread.tlist[i+1]]
                y1 = integrate.odeint(fopdt,ym[i],ts,args=(uf,Km,taum,thetam))
                ym[i+1] = y1[-1]
            return ym

        
        TuningThread.kp = TuningWindow.gainp.value() 
        TuningThread.taup = TuningWindow.taup.value() 
        TuningThread.thetap = TuningWindow.tetap.value()

        if TuningThread.thetap < 0:
            TuningThread.thetap = 0

        x1 = np.zeros(3)
        x1[0] = TuningThread.kp #2.0 # Km
        x1[1] = TuningThread.taup #3.0 # taum
        x1[2] = TuningThread.thetap # thetam
        ym1 = sim_model(x1)

        TuningWindow.dados.setData(TuningThread.tlist , TuningThread.pv)
        TuningWindow.guess.setData(TuningThread.tlist , ym1)
        TuningWindow.fopdt_optimized.clear()


        if TuningThread.thetap <= 0:  #regras de controle imc agressivo // agressive imc control rules
            WorkerThread1.kc = (1/TuningThread.kp)*((TuningThread.taup+0.5*TuningThread.thetap)/(0.1*TuningThread.taup+0.5*TuningThread.thetap))
            WorkerThread1.tauI = TuningThread.taup + 0.5*TuningThread.thetap
            WorkerThread1.tauD = (TuningThread.taup*TuningThread.thetap)/(2*TuningThread.taup+TuningThread.thetap)

        else:   #regras de ziegler-nichols // ziegler-nichols rules
            WorkerThread1.kc = (1.2*TuningThread.taup)/(TuningThread.thetap*TuningThread.kp)
            WorkerThread1.tauI = 2*TuningThread.thetap
            WorkerThread1.tauD = TuningThread.thetap/2

        #atulaizar valores na interface // update values on GUI
        TuningWindow.gainp.setValue(TuningThread.kp)
        TuningWindow.taup.setValue(TuningThread.taup)
        TuningWindow.tetap.setValue(TuningThread.thetap)
        TuningWindow.gainc.setValue(WorkerThread1.kc)
        TuningWindow.tauc.setValue(WorkerThread1.tauI)
        TuningWindow.tetac.setValue(WorkerThread1.tauD)


    def otimizar(self):

        data = TuningThread.pv
        t = TuningThread.tlist
        u = TuningThread.output
        dataexp = [x - TuningThread.pv[0] for x in data]
            
        yp = dataexp
        u0 = u[0]
        yp0 = yp[0]

        # specify number of steps // especificar numero de passos
        ns = len(t)
        delta_t = t[1]-t[0]
        # create linear interpolation of the u data versus time // criar interpolação linear dos dados saida vs t
        uf = interpolate.interp1d(t,u)

        # define first-order plus dead-time approximation // definir aproximação de modelo FOPDT
        def fopdt(y,t,uf,Km,taum,thetam):
            # arguments
            #  y      = output
            #  t      = time
            #  uf     = input linear function (for time shift)
            #  Km     = model gain
            #  taum   = model time constant
            #  thetam = model time constant
            # time-shift u
            try:
                if (t-thetam) <= 0:
                    um = uf(0.0)
                else:
                    um = uf(t-thetam)
            except:
                #print('Error with time extrapolation: ' + str(t))
                um = u0
            # calculate derivative
            dydt = (-(y-yp0) + Km * (um-u0))/taum
            return dydt


        # simulate FOPDT model with x=[Km,taum,thetam]
        def sim_model(x):
            # input arguments
            Km = x[0]
            taum = x[1]
            thetam = x[2]
            # storage for model values
            ym = np.zeros(ns)  # model
            # initial condition
            ym[0] = yp0
            # loop through time steps    
            for i in range(0,ns-1):
                ts = [t[i],t[i+1]]
                y1 = integrate.odeint(fopdt,ym[i],ts,args=(uf,Km,taum,thetam))
                ym[i+1] = y1[-1]
            return ym

        # define objective // definir o objetivo 
        def objective(x):
            # simulate model
            ym = sim_model(x)
            # calculate objective
            obj = 0.0
            for i in range(len(ym)):
                obj = obj + (ym[i]-yp[i])**2    
            # return result
            return obj

        # initial guesses // chutes iniciais
        x0 = np.zeros(3)
        x0[0] = TuningThread.kp #2.0 # Km
        x0[1] = TuningThread.taup #3.0 # taum
        x0[2] = TuningThread.thetap # thetam

        solution = optimize.minimize(objective,x0 , tol = 10)

        # Another way to solve: with bounds on variables
        #bnds = ((0.4, 0.6), (1.0, 10.0), (0.0, 30.0))
        #solution = minimize(objective,x0,bounds=bnds,method='SLSQP')
        x = solution.x

        if x[2] < 0:
            x[2] = 0

        # calculate model with updated parameters // calcular modelo com parametros atualizados
        ym1 = sim_model(x0)
        ym2 = sim_model(x)

        ym1 = [x + TuningThread.pv[0] for x in ym1]
        ym2 = [x + TuningThread.pv[0] for x in ym2] 

        TuningThread.kp = x[0]
        TuningThread.taup = x[1]
        TuningThread.thetap = x[2]

        TuningWindow.dados.setData(t , data)
        TuningWindow.guess.setData(t , ym1)
        TuningWindow.fopdt_optimized.setData(t , ym2)


        if TuningThread.thetap <= 0:  
            WorkerThread1.kc = (1/TuningThread.kp)*((TuningThread.taup+0.5*TuningThread.thetap)/(0.1*TuningThread.taup+0.5*TuningThread.thetap))
            WorkerThread1.tauI = TuningThread.taup + 0.5*TuningThread.thetap
            WorkerThread1.tauD = (TuningThread.taup*TuningThread.thetap)/(2*TuningThread.taup+TuningThread.thetap)

        else:      
            WorkerThread1.kc = (1.2*TuningThread.taup)/(TuningThread.thetap*TuningThread.kp)
            WorkerThread1.tauI = 2*TuningThread.thetap
            WorkerThread1.tauD = TuningThread.thetap/2
       
        TuningWindow.gainp.setValue(TuningThread.kp)
        TuningWindow.taup.setValue(TuningThread.taup)
        TuningWindow.tetap.setValue(TuningThread.thetap)
        TuningWindow.gainc.setValue(WorkerThread1.kc)
        TuningWindow.tauc.setValue(WorkerThread1.tauI)
        TuningWindow.tetac.setValue(WorkerThread1.tauD)
        TuningWindow.squad_init.setValue(objective(x0))
        TuningWindow.squad_final.setValue(objective(x))

class WorkerThread1(QThread):  # Thread do controlador // Controller thread
    
    #Dev2 = sistema simulado 
    #Dev3 = bancada
    dev = "Dev2"  # device of communication  // dispositivo de comunicação
    input_h = "ai3" # input channel
    output = "ao0" # output channel
    var=1
    var_signal_filter = False
    var_derivative_filter = False
    tempo = []
    change_value = pyqtSignal(int)
    change_value2 = pyqtSignal(list)
    malha_aberta = pyqtSignal(int)
    malha_fechada = pyqtSignal(int)
    sp = []
    nivel = []
    bomba = []
    indice = 0
    kc = 12.49
    tauI = 5.76
    tauD = 1.44
    wr = 0
    dt = 0
    te = 0
    IAE = 0
    ITAE = 0
    ISE = 0
    IE = 0
    ie = []  # integral of the error
    filter_list = []
    smoothfilter_list = []

    def run(self):
        #abrir workbook do excel, e criar uma folha // open excel and creat a sheet
        wb = Workbook()
        sheet1 = wb.add_sheet('Sheet 1')
        data_e_hora_atuais = datetime.now()
        data_e_hora = data_e_hora_atuais.strftime('%d-%m-%Y_%H_%M')
        #local onde os dados serão salvos // path to save data
        local_save= ('C:\\Users\\luans\\Desktop\\'+data_e_hora+'.xls')

        sheet1.write(1, 0, 'tempo(s)')
        sheet1.write(1, 1, 'nível(volts)')
        sheet1.write(1, 2, 'nível(mm)')
        sheet1.write(1, 3, 'bomba(%)')
        sheet1.write(1, 4, 'bomba(volts)')
        sheet1.write(1, 5, 'hora')
        sheet1.write(1, 6, 'setpoint(mm)')
        sheet1.write(1, 7, 'IE')
        sheet1.write(1, 8, 'IAE')
        sheet1.write(1, 9, 'ISE')
        sheet1.write(1, 10, 'ITAE')

        e = []   # error
        dpv = [] # derivative of the pv
        P = []   # proportional
        I = []   # integral
        D = []   # derivative
        xi = []
        dxi = []
        fd = []
        PID = 0
        hf = 0

        t=0
        dt=1
        WorkerThread1.dt = dt
        
        
        with nidaqmx.Task() as task, nidaqmx.Task() as task1:
            task.ai_channels.add_ai_voltage_chan(WorkerThread1.dev+"/"+WorkerThread1.input_h, terminal_config = TerminalConfiguration.RSE)
            task1.ao_channels.add_ao_voltage_chan(WorkerThread1.dev+"/"+WorkerThread1.output)
            y=[]
            z=[]
            i=2
            c=0
            vf=True
            while vf == True:
                if WorkerThread1.var == 4: # PID Control
                    #ler dados // read data
                    data = task.read(number_of_samples_per_channel=1)
                    h = CalibrationWindow.a*data[0]+CalibrationWindow.b
                    if h<0:
                        h = 0

                    #filtro de sinal // signal filter
                    if WorkerThread1.var_signal_filter == True:  
                        if len(WorkerThread1.filter_list)>=10:
                            WorkerThread1.filter_list.append(h)
                            del  WorkerThread1.filter_list[0] 
                            WorkerThread1.smoothfilter_list = signal.savgol_filter(WorkerThread1.filter_list, 7, 2)
                            hf = WorkerThread1.smoothfilter_list[4]
                                    
                        else:
                            WorkerThread1.filter_list.append(h)
                            hf = h

                    elif WorkerThread1.var_signal_filter == False:
                        if len(WorkerThread1.filter_list)>=5:
                            WorkerThread1.filter_list.append(h)
                            del  WorkerThread1.filter_list[0]
                            hf = h 
                                    
                        else:
                            WorkerThread1.filter_list.append(h)
                            hf = h

                    if hf<0:
                        hf = 0

                    WorkerThread1.sp.append(MainWindow.setpoint)  
                    WorkerThread1.nivel.append(hf)

                    e.append(WorkerThread1.sp[c] - hf)
    
                    if c >= 1:  # calculo começando na segunda iteração // calculation starts in the second iteration
                        dpv.append((WorkerThread1.nivel[c]-WorkerThread1.nivel[c-1])/dt)
                        
                        WorkerThread1.ie.append(WorkerThread1.ie[c-1] + e[c] * dt)  
                        #calcular indices de desemepenho // calculate performance indices
                        WorkerThread1.IE = WorkerThread1.ie[c] # Error Integral
                        WorkerThread1.IAE = WorkerThread1.IAE + abs(e[c]) * dt # Integral of the absolute error
                        WorkerThread1.ISE = WorkerThread1.ISE + e[c]**2 * dt # Integral of the squared error
                        WorkerThread1.ITAE = WorkerThread1.ITAE + abs(e[c]) * WorkerThread1.te * dt # Integral of the  time absolute error
                        WorkerThread1.te = WorkerThread1.te + dt
                    else:
                        dpv.append(0)
                        WorkerThread1.ie.append(0)
                    
                    # resposta do PID // PID output   
                    P.append(WorkerThread1.kc * e[c]) 
                    I.append(WorkerThread1.kc/WorkerThread1.tauI * WorkerThread1.ie[c])  
                    D.append( -WorkerThread1.kc * WorkerThread1.tauD * dpv[c])
                    PID = xi[0] + P[c] + I[c] + D[c] 
                    dxi.append((PID-xi[c-1])/dt)
                    fd.append(TuningThread.alfa * WorkerThread1.tauD * dxi[c])

                    if c >= 1: 
                        if WorkerThread1.var_derivative_filter == True:
                            xi.append(xi[0] + P[c] + I[c] + D[c] - fd[c])
                        elif WorkerThread1.var_derivative_filter == False:
                            xi.append(xi[0] + P[c] + I[c] + D[c])

                    else:
                        dxi.append(0)
                        fd.append(0)
                        xi.append(0)

                    if xi[c] > 90:
                        WorkerThread1.wr = 90
                        I[WorkerThread1.indice+1] = I[WorkerThread1.indice] #integral reset
                    elif xi[c] < 0:
                        WorkerThread1.wr = 0
                        I[WorkerThread1.indice+1] = I[WorkerThread1.indice] #integral reset
                    else:
                        WorkerThread1.wr = xi[c]

                    WorkerThread1.bomba.append(WorkerThread1.wr)
                    task1.write(WorkerThread1.wr/20)

                    #obter hora do momento da aquisição da amostra // time of data acquisition
                    hora=datetime.now()
                    tempo_decorrido=hora.strftime('%H:%M:%S')
                    WorkerThread1.tempo.append(t)
                    
                    #escrever dados na folha criada no workbook do excel // write data on excel sheet
                    sheet1.write(i, 0, t) #(linha, coluna, dados) // (row, column, data)
                    sheet1.write(i, 1, data[0])
                    sheet1.write(i, 2, hf)
                    sheet1.write(i, 3, WorkerThread1.wr)
                    sheet1.write(i, 4, WorkerThread1.wr/20)
                    sheet1.write(i, 5, tempo_decorrido)
                    sheet1.write(i, 6, WorkerThread1.sp[c])
                    sheet1.write(i, 7, WorkerThread1.IE)
                    sheet1.write(i, 8, WorkerThread1.IAE)
                    sheet1.write(i, 9, WorkerThread1.ISE)
                    sheet1.write(i, 10, WorkerThread1.ITAE)

                    WorkerThread1.indice = c
                    t=t+dt
                    i=i+1
                    c=c+1
                    
                    z.append(WorkerThread1.wr)
                    y.append(hf)
                    
                    self.change_value.emit(y)
                    self.change_value2.emit(z)
                    self.malha_fechada.emit(1)
                    
                    wb.save(local_save)
                    time.sleep(dt)

                elif WorkerThread1.var == 1:
                    #ler dados // read data
                    data = task.read(number_of_samples_per_channel=1)
                    h = CalibrationWindow.a*data[0]+CalibrationWindow.b
                    if h<0:
                        h = 0
                    
                    #filtro de sinal // signal filter
                    if WorkerThread1.var_signal_filter == True:  
                        if len(WorkerThread1.filter_list)>=5:
                            WorkerThread1.filter_list.append(h)
                            del  WorkerThread1.filter_list[0] 
                            WorkerThread1.smoothfilter_list = signal.savgol_filter(WorkerThread1.filter_list, 5, 2)
                            hf = WorkerThread1.smoothfilter_list[4]
                                    
                        else:
                            WorkerThread1.filter_list.append(h)
                            hf = h 
                    elif WorkerThread1.var_signal_filter == False:
                        if len(WorkerThread1.filter_list)>=5:
                            WorkerThread1.filter_list.append(h)
                            del  WorkerThread1.filter_list[0]
                            hf = h 
                                    
                        else:
                            WorkerThread1.filter_list.append(h)
                            hf = h
                    
                    if hf<0:
                        hf = 0

                    WorkerThread1.sp.append(MainWindow.setpoint)  
                    WorkerThread1.nivel.append(hf)

                    e.append(hf - hf)
    
                    if c >= 1:  # calculo começando na segunda iteração
                        dpv.append((WorkerThread1.nivel[c]-WorkerThread1.nivel[c-1])/dt) 
                        WorkerThread1.ie.append(WorkerThread1.ie[c-1] + e[c] * dt)  
                    else:
                        dpv.append(0)
                        WorkerThread1.ie.append(0)
                    P.append(WorkerThread1.kc * e[c]) 
                    I.append(WorkerThread1.kc/WorkerThread1.tauI * WorkerThread1.ie[c])  
                    D.append( -WorkerThread1.kc * WorkerThread1.tauD * dpv[c])

                    dxi.append(0)
                    fd.append(TuningThread.alfa * WorkerThread1.tauD * dxi[c])

                    if c >= 1:
                        xi.append(xi[0] + P[c] + I[c] + D[c])
    
                    else:
                        xi.append(0)

                    WorkerThread1.bomba.append(WorkerThread1.wr)
                    task1.write(WorkerThread1.wr/20)
                    if TuningThread.getD1 == True:
                        TuningThread.t1 = TuningThread.count
                        TuningThread.getD1 = False
                    #obter hora do momento da aquisição da amostra
                    hora=datetime.now()
                    tempo_decorrido=hora.strftime('%H:%M:%S')
                    WorkerThread1.tempo.append(t)
                    #escrever dados na folha criada no workbook do excel
                    sheet1.write(i, 0, t) #(linha, coluna, dados)
                    sheet1.write(i, 1, data[0])
                    sheet1.write(i, 2, hf)
                    sheet1.write(i, 3, WorkerThread1.wr)
                    sheet1.write(i, 4, WorkerThread1.wr/20)
                    sheet1.write(i, 5, tempo_decorrido)
                    sheet1.write(i, 6, WorkerThread1.sp[c])
                    WorkerThread1.indice = c
                    t=t+dt
                    i=i+1
                    c=c+1

                    z.append(WorkerThread1.wr)
                    y.append(hf)
                    
                    self.malha_aberta.emit(1)
                    
                    wb.save(local_save)
                    time.sleep(dt)


                elif WorkerThread1.var == 2:
                    task1.write(0)
                    time.sleep(dt)

                elif WorkerThread1 == 3:
                    vf=False
          
##########################################################################################################        

## Thread principal // Main thread (GUI thread)
class MainWindow(QtWidgets.QMainWindow): 
    switch_window = QtCore.pyqtSignal()
    open_window = QtCore.pyqtSignal()
    setpoint= 0
    btn = 0

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.title = 'Bancada de controle de nível'
        self.left = 100
        self.top = 100
        self.width = 1200
        self.height = 600
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.home()

    def closeEvent(self, event):
        with nidaqmx.Task() as task2:
            task2.ao_channels.add_ao_voltage_chan(WorkerThread1.dev+"/"+WorkerThread1.output)
            task2.write(0) 
            
    def home(self):
        self.showMaximized()
        layout_principal = QtWidgets.QHBoxLayout()
        dual_left = QtWidgets.QHBoxLayout()
        dual_right = QtWidgets.QHBoxLayout()
        descricao = QtWidgets.QVBoxLayout()
        control_panel = QtWidgets.QVBoxLayout()

        vbox = QtWidgets.QVBoxLayout()
        hbox = QtWidgets.QHBoxLayout()
        self.graphWidget = pg.PlotWidget()  # nivel
        self.graphWidget.setLabel('left', "Nivel (mm)")
        self.graphWidget.setLabel('bottom', "Tempo (s)")
        self.graphWidget2 = pg.PlotWidget()  # bomba
        self.graphWidget2.setLabel('left', "Vazão da bomba (%)")
        self.graphWidget2.setLabel('bottom', "Tempo (s)")
        self.btn = QtWidgets.QPushButton("INICIAR")
        self.btn2 = QtWidgets.QPushButton("SINTONIA")
        self.btn3 = QtWidgets.QPushButton("PAUSAR")
        self.btn4 = QtWidgets.QPushButton("RESUMIR")
        central_widget = QtWidgets.QWidget()
        central_widget.setLayout(layout_principal)
        self.btn.clicked.connect(self.startPlotAq)
        self.btn2.clicked.connect(self.tuning)
        self.btn3.clicked.connect(self.pausePlot)
        self.btn4.clicked.connect(self.resumePlot)
        #self.btn.clicked.connect(self.startPlotWr)
     
        a1 = QtWidgets.QHBoxLayout()
        self.labelsp = QtWidgets.QLabel("Set-point:", self)
        self.out = QtWidgets.QDoubleSpinBox()
        self.out.setRange(0,120)
        a1.addWidget(self.labelsp)
        a1.addWidget(self.out)

        a2 = QtWidgets.QHBoxLayout()
        self.labelbomba = QtWidgets.QLabel("Vazão da Bomba (%):", self)
        self.bomba = QtWidgets.QDoubleSpinBox()
        self.bomba.setRange(0,200)
        a2.addWidget(self.labelbomba)
        a2.addWidget(self.bomba)

        a7 = QtWidgets.QHBoxLayout()
        self.labelnivel = QtWidgets.QLabel("Nível (mm):", self)
        self.nivelspin = QtWidgets.QDoubleSpinBox()
        self.nivelspin.setRange(-100000,1000000)
        a7.addWidget(self.labelnivel)
        a7.addWidget(self.nivelspin)

        a3 = QtWidgets.QHBoxLayout()
        self.labelgain = QtWidgets.QLabel("Ganho do Controlador (Kc):", self)
        self.gainspin = QtWidgets.QDoubleSpinBox()
        self.gainspin.setRange(0,100000)
        a3.addWidget(self.labelgain)
        a3.addWidget(self.gainspin)

        a4 = QtWidgets.QHBoxLayout()
        self.labeltint = QtWidgets.QLabel("Tempo Integral (Ti):", self)
        self.tintspin = QtWidgets.QDoubleSpinBox()
        self.tintspin.setRange(0.01,100000)
        a4.addWidget(self.labeltint)
        a4.addWidget(self.tintspin)

        a5 = QtWidgets.QHBoxLayout()
        self.labeltd = QtWidgets.QLabel("Tempo Derivativo (Td):", self)
        self.tdspin = QtWidgets.QDoubleSpinBox()
        self.tintspin.setRange(0.0,100000)
        a5.addWidget(self.labeltd)
        a5.addWidget(self.tdspin)

        a6 = QtWidgets.QHBoxLayout()
        self.pidcheck = QtWidgets.QCheckBox("Controle PID")
        a6.addWidget(self.pidcheck)

        performance = QtWidgets.QVBoxLayout()

        a12 = QtWidgets.QHBoxLayout()
        self.labelper = QtWidgets.QLabel("Indices de Performance:", self)
        a12.addWidget(self.labelper)

        a8 = QtWidgets.QHBoxLayout()
        self.labelie = QtWidgets.QLabel("IE:", self)
        self.iespin = QtWidgets.QDoubleSpinBox()
        self.iespin.setRange(-10000000,10000000)
        self.iespin.setDecimals(5)
        a8.addWidget(self.labelie)
        a8.addWidget(self.iespin)

        a9 = QtWidgets.QHBoxLayout()
        self.labeliae = QtWidgets.QLabel("IAE:", self)
        self.iaespin = QtWidgets.QDoubleSpinBox()
        self.iaespin.setRange(-10000000,10000000)
        self.iaespin.setDecimals(5)
        a9.addWidget(self.labeliae)
        a9.addWidget(self.iaespin)

        a10 = QtWidgets.QHBoxLayout()
        self.labelise = QtWidgets.QLabel("ISE:", self)
        self.isespin = QtWidgets.QDoubleSpinBox()
        self.isespin.setRange(-10000000,10000000)
        self.isespin.setDecimals(5)
        a10.addWidget(self.labelise)
        a10.addWidget(self.isespin)

        a11 = QtWidgets.QHBoxLayout()
        self.labelitae = QtWidgets.QLabel("ITAE:", self)
        self.itaespin = QtWidgets.QDoubleSpinBox()
        self.itaespin.setRange(-10000000,10000000)
        self.itaespin.setDecimals(5)
        a11.addWidget(self.labelitae)
        a11.addWidget(self.itaespin)

        a13 = QtWidgets.QHBoxLayout()
        self.labelb = QtWidgets.QLabel("Instruções: \nO experimento será iniciado em malha aberta ao clicar em 'INICIAR'.\nCaso esteja em malha aberta a vazão da bomba pode\nser alterada no seu respectivo campo de preenchimento", self)
        self.labelb.setFont(QtGui.QFont('Arial', 10))
        a13.addWidget(self.labelb)

        a14 = QtWidgets.QHBoxLayout()
        self.labelb1 = QtWidgets.QLabel("Ao marcar o campo 'Controle PID', a malha será fechada utilizando\no controle do tipo PID.\nNesse modo o campo 'Set-point', e os parametros 'Kc', 'Ti' e 'Td'\npodem ser alterados livremente.\nCaso queira ajustar os parâmetros utilizando o método de\nZiegler-Nichols clique em 'SINTONIA'", self)
        self.labelb1.setFont(QtGui.QFont('Arial', 10))
        a14.addWidget(self.labelb1)

        a15 = QtWidgets.QHBoxLayout()
        self.labelb2 = QtWidgets.QLabel(" ", self)
        a15.addWidget(self.labelb2)

        a16 = QtWidgets.QHBoxLayout()
        self.signal_filter = QtWidgets.QCheckBox("Filtro de sinal")
        a16.addWidget(self.signal_filter)

        a17 = QtWidgets.QHBoxLayout()
        self.derivative_filter = QtWidgets.QCheckBox("Filtro derivativo")
        a17.addWidget(self.derivative_filter)

        #control_panel
        control_panel.addLayout(a7)
        control_panel.addLayout(a1)
        control_panel.addLayout(a2)
        control_panel.addLayout(a3)
        control_panel.addLayout(a4)
        control_panel.addLayout(a5)
        control_panel.addLayout(a6)
        control_panel.addLayout(a16)
        control_panel.addLayout(a17)
        control_panel.addLayout(a12)
        control_panel.addLayout(a8)
        control_panel.addLayout(a9)
        control_panel.addLayout(a10)
        control_panel.addLayout(a11)
        control_panel.addLayout(a13)
        control_panel.addLayout(a14)
        control_panel.addLayout(a15)

        #dual_right
        dual_right.addLayout(control_panel)

        #layout_principal
        layout_principal.addLayout(dual_left)
        layout_principal.addLayout(dual_right)

        #dual_left
        vbox.addWidget(self.graphWidget)
        vbox.addWidget(self.graphWidget2)
        hbox.addWidget(self.btn)
        hbox.addWidget(self.btn2)
        #hbox.addWidget(self.btn3)      #####
        #hbox.addWidget(self.btn4)      #####
        vbox.addLayout(hbox)
        dual_left.addLayout(vbox)
        self.setCentralWidget(central_widget)
        self.show()
        
        self.pidcheck.clicked.connect(self.pid_check)
        self.derivative_filter.clicked.connect(self.derivative_filter_check)
        self.signal_filter.clicked.connect(self.signal_filter_check)
        self.out.editingFinished.connect(self.on_value_changed)
        self.bomba.editingFinished.connect(self.bomba_value_changed)
        self.tintspin.editingFinished.connect(self.tint_value_changed)
        self.tdspin.editingFinished.connect(self.td_value_changed)
        self.gainspin.editingFinished.connect(self.gain_value_changed)

        self.nivel = pg.PlotCurveItem(clear=True, pen = pg.mkPen('r', width=2) , name=  'Nível')
        self.setpointgraph = pg.PlotCurveItem(clear=True, pen = pg.mkPen('b', width=2) , name='Setpoint')
        self.graphWidget.addLegend()

        self.graphWidget.addItem(self.nivel)
        self.graphWidget.addItem(self.setpointgraph)
        
        self.pidcheck.setEnabled(False)
        self.derivative_filter.setEnabled(False)

             
    def pid_check(self): # turn PID control on/off
        if self.pidcheck.isChecked() == True:
            WorkerThread1.ie[WorkerThread1.indice] = 0
            WorkerThread1.var = 4
            MainWindow.setpoint = WorkerThread1.nivel[WorkerThread1.indice]
            self.out.setValue(WorkerThread1.nivel[WorkerThread1.indice])
            self.btn2.setEnabled(False)
        elif self.pidcheck.isChecked() == False:
            WorkerThread1.var = 1
            self.btn2.setEnabled(True)

    def signal_filter_check(self): # turn signal filter on/off
        if self.signal_filter.isChecked() == True:
            WorkerThread1.var_signal_filter = True
        elif self.pidcheck.isChecked() == False:
            WorkerThread1.var_signal_filter = False

    def derivative_filter_check(self): #turn derivative filter on/off
        if self.derivative_filter.isChecked() == True:
            WorkerThread1.var_derivative_filter = True
            print("True")
        elif self.derivative_filter.isChecked() == False:
            WorkerThread1.var_derivative_filter = False
            print("False")

    def tuning(self):
        self.open_window.emit()

    def resumePlot(self):
        WorkerThread1.var = 1  
        WorkerThread1.dt=1

    def pausePlot(self):
        WorkerThread1.var = 2  
        

    def on_value_changed(self):    
        MainWindow.setpoint = self.out.value()
        WorkerThread1.te = 0
        WorkerThread1.ie[WorkerThread1.indice] = 0
        WorkerThread1.IAE = 0
        WorkerThread1.ISE = 0
        WorkerThread1.ITAE = 0

    def bomba_value_changed(self):
        WorkerThread1.wr = self.bomba.value()
           
    def tint_value_changed(self):
        WorkerThread1.tauI = self.tintspin.value()

    def td_value_changed(self):
        WorkerThread1.tauD = self.tdspin.value()

    def gain_value_changed(self):
        WorkerThread1.kc = self.gainspin.value()       

    def startPlotAq(self):  #aquisição de dados pela worker thread // call worker thread
        self.thread = WorkerThread1()
        self.thread2 = TuningThread()
        
        #malha_fechada // closed loop
        self.thread.malha_fechada.connect(self.setValueNivel)
        self.thread.malha_fechada.connect(self.setValueSetpoint)
        self.thread.malha_fechada.connect(self.setValueBomba)
        self.thread.malha_fechada.connect(self.setInfoFechada)

        #malha_aberta // open loop
        self.thread.malha_aberta.connect(self.setValueNivel)
        self.thread.malha_aberta.connect(self.setValueBomba)
        self.thread.malha_aberta.connect(self.setInfoAberta)

        self.thread.start()
        self.btn.setEnabled(False)
        self.pidcheck.setEnabled(True)
        self.derivative_filter.setEnabled(True)

    def setInfoAberta(self, dados):
        self.setpointgraph.setData([0],[0])
        self.gainspin.setValue(WorkerThread1.kc)
        self.tintspin.setValue(WorkerThread1.tauI)
        self.tdspin.setValue(WorkerThread1.tauD)


    def setInfoFechada(self, dados):
        self.bomba.setValue(WorkerThread1.bomba[WorkerThread1.indice])
        self.iespin.setValue(WorkerThread1.IE)
        self.isespin.setValue(WorkerThread1.ISE)
        self.iaespin.setValue(WorkerThread1.IAE)
        self.itaespin.setValue(WorkerThread1.ITAE)
    
    def setValueSetpoint(self, nivel):
        self.setpointgraph.setData(WorkerThread1.tempo, WorkerThread1.sp)

    def setValueNivel(self, nivel):
        self.nivel.setData(WorkerThread1.tempo, WorkerThread1.nivel)
        self.nivelspin.setValue(WorkerThread1.nivel[WorkerThread1.indice])
    
    def setValueBomba(self, bomba):
        self.graphWidget2.plot(WorkerThread1.tempo,WorkerThread1.bomba, clear=True, pen = pg.mkPen('w', width=2))


###################################################################################################   

## Interface de sintonia do controlador // Controller tuning GUI

class TuningWindow(QtWidgets.QMainWindow):
    switch_window = QtCore.pyqtSignal()
    data1 = QtCore.pyqtSignal(int)

    graph1 = 0
    dados = 0
    guess = 0
    fopdt_optimized = 0
    graph2 = 0
    gainp = 0
    gainc = 0
    tauc = 0
    taup = 0
    tetac = 0
    tetap = 0
    plot = 0
    identification_type = ""
    squad_init = 0
    squad_final = 0

    def closeEvent(self, event):
        self.tuningthread = TuningThread()
        self.tuningthread.exit()

    def __init__(self, *args, **kwargs):
        super(TuningWindow, self).__init__(*args, **kwargs)

        self.title = 'Sintonia de Parâmetros'
        self.left = 100
        self.top = 100
        self.width = 1200
        self.height = 600
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.home()


    def home(self):
        self.showMaximized()
        vbox = QtWidgets.QVBoxLayout()
        hbox = QtWidgets.QHBoxLayout()
        dual = QtWidgets.QHBoxLayout()
        dual_vert = QtWidgets.QVBoxLayout()
        self.graphWidget = pg.PlotWidget()  # nivel
        self.graphWidget.setLabel('left', "Nivel (mm)")
        self.graphWidget.setLabel('bottom', "Tempo (s)")
        self.graphWidget2 = pg.PlotWidget()  # bomba
        self.graphWidget2.setLabel('left', "Vazão da bomba (%)")
        self.graphWidget2.setLabel('bottom', "Tempo (s)")
        self.btn = QtWidgets.QPushButton("INICIAR")
        self.btn2 = QtWidgets.QPushButton("DEGRAU +10%")
        self.btn3 = QtWidgets.QPushButton("DEGRAU -10%")
        self.btn4 = QtWidgets.QPushButton("FINALIZAR")
        self.btn6 = QtWidgets.QPushButton("SIMULAR MODELO")
        self.btn7 = QtWidgets.QPushButton("MÉTODO DA TANGENTE")
        self.btn5 = QtWidgets.QPushButton("OTIMIZAÇÃO")
        central_widget = QtWidgets.QWidget()
        central_widget.setLayout(vbox)
        self.btn.clicked.connect(self.startPlotAq)
        self.btn2.clicked.connect(self.degrauPlus)
        self.btn3.clicked.connect(self.degrauMinus)
        self.btn4.clicked.connect(self.stop)
        #self.btn.clicked.connect(self.startPlotWr)
        dual_vert.addWidget(self.graphWidget)
        dual_vert.addWidget(self.graphWidget2)
        dual.addLayout(dual_vert)
        #vbox.addWidget(self.out)
        hbox.addWidget(self.btn)
        hbox.addWidget(self.btn2)
        hbox.addWidget(self.btn3)
        hbox.addWidget(self.btn4)
        hbox.addWidget(self.btn6)
        hbox.addWidget(self.btn7)
        hbox.addWidget(self.btn5)
        
        '''
        ident_type = QtWidgets.QHBoxLayout()
        self.labelidentification = QtWidgets.QLabel("Método de identificação:", self)
        self.identification_combo = QtWidgets.QComboBox()
        self.identification_combo.addItems(["Método da Tangente","Otimização"])
        ident_type.addWidget(self.labelidentification)
        ident_type.addWidget(self.identification_combo)
        '''
        
        processo = QtWidgets.QVBoxLayout()
        gain = QtWidgets.QHBoxLayout()
        self.labelgain = QtWidgets.QLabel("Ganho:", self) 
        self.spingain = QtWidgets.QDoubleSpinBox()
        self.spingain.setRange(-10000000,10000000)
        TuningWindow.gainp = self.spingain
        gain.addWidget(self.labelgain)
        gain.addWidget(self.spingain)
        
        tconst = QtWidgets.QHBoxLayout()
        self.labeltconst = QtWidgets.QLabel("Constante de Tempo:", self)  
        self.spintconst = QtWidgets.QDoubleSpinBox()  
        TuningWindow.taup = self.spintconst
        self.spintconst.setRange(-10000000,10000000)
        tconst.addWidget(self.labeltconst)
        tconst.addWidget(self.spintconst)

        tdelay = QtWidgets.QHBoxLayout()
        self.labeldelay = QtWidgets.QLabel("Tempo Morto:", self)   
        self.spindelay = QtWidgets.QDoubleSpinBox()
        TuningWindow.tetap = self.spindelay
        self.spindelay.setRange(-10000000,10000000)
        tdelay.addWidget(self.labeldelay)
        tdelay.addWidget(self.spindelay)

        squadrado_init = QtWidgets.QHBoxLayout()
        self.squadrado1l = QtWidgets.QLabel("ΣE\u00b2 inicial:", self)   
        self.squadrado1spin = QtWidgets.QDoubleSpinBox()
        TuningWindow.squad_init = self.squadrado1spin
        self.squadrado1spin.setRange(-10000000,10000000)
        squadrado_init.addWidget(self.squadrado1l)
        squadrado_init.addWidget(self.squadrado1spin)

        squadrado_final = QtWidgets.QHBoxLayout()
        self.squadrado2l = QtWidgets.QLabel("ΣE\u00b2 final:", self)   
        self.squadrado2spin = QtWidgets.QDoubleSpinBox()
        TuningWindow.squad_final = self.squadrado2spin
        self.squadrado2spin.setRange(-10000000,10000000)
        squadrado_final.addWidget(self.squadrado2l)
        squadrado_final.addWidget(self.squadrado2spin)

        a3 = QtWidgets.QHBoxLayout()
        self.labelgain = QtWidgets.QLabel("Ganho do Controlador(Kc):", self)
        self.gainspin = QtWidgets.QDoubleSpinBox()
        TuningWindow.gainc = self.gainspin
        self.gainspin.setRange(-10000000,10000000)
        a3.addWidget(self.labelgain)
        a3.addWidget(self.gainspin)

        a4 = QtWidgets.QHBoxLayout()
        self.labeltint = QtWidgets.QLabel("Tempo Integral(Ti):", self)
        self.tintspin = QtWidgets.QDoubleSpinBox()
        TuningWindow.tauc = self.tintspin
        self.tintspin.setRange(-10000000,10000000)
        a4.addWidget(self.labeltint)
        a4.addWidget(self.tintspin)

        a5 = QtWidgets.QHBoxLayout()
        self.labeltd = QtWidgets.QLabel("Tempo Derivativo(Td):", self)
        self.tdspin = QtWidgets.QDoubleSpinBox()
        self.tdspin.setRange(-10000000,10000000)
        TuningWindow.tetac = self.tdspin
        a5.addWidget(self.labeltd)
        a5.addWidget(self.tdspin)

        a6 = QtWidgets.QHBoxLayout()
        self.instr = QtWidgets.QLabel("Instruções: \nPara iniciar o experimento de sintonia com teste degrau clique em 'INICIAR'. Escolha se a entrada degrau da vazão da bomba será positiva ou negativa clicando em 'DEGRAU +10%' ou 'DEGRAU -10%'\nQuando o processo atingir o estado estacionário, clique em 'FINALIZAR' para encerrar o experimento e calcular os parâmetros do modelo do processo e do controlador.", self)
        self.instr.setFont(QtGui.QFont('Arial', 10))
        a6.addWidget(self.instr)

        #processo.addLayout(ident_type)
        processo.addLayout(gain)
        processo.addLayout(tconst)
        processo.addLayout(tdelay)
        processo.addLayout(squadrado_init)
        processo.addLayout(squadrado_final)
        processo.addLayout(a3)
        processo.addLayout(a4)
        processo.addLayout(a5)
        vbox.addLayout(a6)


        dual.addLayout(processo)
        vbox.addLayout(dual)
        vbox.addLayout(hbox)
        self.setCentralWidget(central_widget)
        self.show()
        
        self.nivel = pg.PlotCurveItem(clear=True, pen="r", name = "Nível")
        self.setpointgraph = pg.PlotCurveItem(clear=True, pen="b", name = "Setpoint")
        

        self.graphWidget.addItem(self.nivel)
        TuningWindow.plot = self.graphWidget
        
        TuningWindow.graph1 = self.graphWidget
        TuningWindow.graph2 = self.graphWidget2 

        TuningWindow.dados = pg.PlotCurveItem(clear=True, pen = pg.mkPen('r', width=2) , name=  'Dados de Processo')
        TuningWindow.guess = pg.PlotCurveItem(clear=True, pen = pg.mkPen('b', width=2) , name='Modelo Simulado')
        TuningWindow.fopdt_optimized = pg.PlotCurveItem(clear=True, pen = pg.mkPen('g', width=2) , name='Otimização')
        TuningWindow.graph1.addLegend()

        TuningWindow.graph1.addItem(TuningWindow.dados)
        TuningWindow.graph1.addItem(TuningWindow.guess)
        TuningWindow.graph1.addItem(TuningWindow.fopdt_optimized)


        self.btn2.clicked.connect(TuningThread.getData1)
        self.btn3.clicked.connect(TuningThread.getData1)
        self.btn6.clicked.connect(TuningThread.update_model)
        self.btn7.clicked.connect(TuningThread.getData2)
        self.btn5.clicked.connect(TuningThread.otimizar)

    def optimization(self):
        self.btn5.setEnabled(False)
        
    def ident(self):
        TuningWindow.identification_type = self.identification_combo.currentText()

    def degrauPlus(self):
        WorkerThread1.wr = WorkerThread1.wr + 10
        TuningThread.Adegrau = 10 
        self.btn2.setEnabled(False)
        self.btn3.setEnabled(False)

    def stop(self):
        TuningThread.chave = False
        
    def set_data(self):
        self.spingain.setValue(TuningThread.kp)
        self.spintconst.setValue(TuningThread.taup)
        self.spindelay.setValue(TuningThread.thetap)
        self.gainspin.setValue(WorkerThread1.kc)
        self.tintspin.setValue(WorkerThread1.tauI)
        self.tdspin.setValue(WorkerThread1.tauD)

    def degrauMinus(self):
        WorkerThread1.wr = WorkerThread1.wr - 10
        TuningThread.Adegrau = -10
        self.btn2.setEnabled(False)
        self.btn3.setEnabled(False)

    def on_value_changed(self):
        
        MainWindow.setpoint = self.out.value()
           

    def startPlotAq(self):  #aquisição de dados pela worker thread // starts data acquisition
        self.thread = TuningThread()
        self.thread.change_value.connect(self.setValue)
        self.thread.change_value.connect(self.setValue2)
        self.thread.start()
        self.btn.setEnabled(False)
        TuningThread.chave = True
        TuningThread.var = 1
        
    def setValue(self, nivel):
        self.nivel.setData(TuningThread.tlist, TuningThread.pv, pen = pg.mkPen('r', width=2))
    
    def setValue2(self, bomba):
        self.graphWidget2.plot(TuningThread.tlist, TuningThread.output, clear=True, pen = pg.mkPen('w', width=2))

    
#################################################################################################3

# Interface de calibração do sensor de nível // Level sensor calibration GUI
class CalibrationWindow(QtWidgets.QWidget):
    switch_window = QtCore.pyqtSignal()
    a=0
    b=0

    def __init__(self, *args, **kwargs):
        super(CalibrationWindow, self).__init__(*args, **kwargs)

        self.title = 'Calibração de sensor de nível'
        self.left = 300
        self.top = 100
        self.width = 300
        self.height = 500
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.home()


    def home(self):
        
        vbox = QtWidgets.QVBoxLayout
        hbox1 = QtWidgets.QHBoxLayout()
        hbox2 = QtWidgets.QHBoxLayout()
        hbox3 = QtWidgets.QHBoxLayout()
        hbox4 = QtWidgets.QHBoxLayout()

        self.label1 = QtWidgets.QLabel("Ponto 1:", self)
        self.label2 = QtWidgets.QLabel("Ponto 2:", self)
        self.label3 = QtWidgets.QLabel("Ponto 3:", self)
        self.label4 = QtWidgets.QLabel("Ponto 4:", self)
        
        alt1=QtWidgets.QLabel("                    Altura [mm]:", self)
        self.p1 = QtWidgets.QDoubleSpinBox()
        self.p1.setRange(-10000,10000)
        volt1=QtWidgets.QLabel("                    Tensão [V]:", self)
        self.p11 = QtWidgets.QDoubleSpinBox()
        self.p11.setRange(-10,10)
        self.p11.setDecimals(4)
        btn1= QtWidgets.QPushButton("LER TENSÃO")

        alt2=QtWidgets.QLabel("                    Altura [mm]:", self)
        self.p2 = QtWidgets.QDoubleSpinBox()
        self.p2.setRange(-10000,10000)
        volt2=QtWidgets.QLabel("                    Tensão [V]:", self)
        self.p22 = QtWidgets.QDoubleSpinBox()
        self.p22.setRange(-10,10)
        self.p22.setDecimals(4)
        btn2= QtWidgets.QPushButton("LER TENSÃO")

        alt3=QtWidgets.QLabel("                    Altura [mm]:", self)
        self.p3 = QtWidgets.QDoubleSpinBox()
        self.p3.setRange(-10000,10000)
        volt3=QtWidgets.QLabel("                    Tensão [V]:", self)
        self.p33 = QtWidgets.QDoubleSpinBox()
        self.p33.setRange(-10,10)
        self.p33.setDecimals(4)
        btn3= QtWidgets.QPushButton("LER TENSÃO")

        alt4=QtWidgets.QLabel("                    Altura [mm]:", self)
        self.p4 = QtWidgets.QDoubleSpinBox()
        self.p4.setRange(-10000,10000)
        volt4=QtWidgets.QLabel("                    Tensão [V]:", self)
        self.p44 = QtWidgets.QDoubleSpinBox()
        self.p44.setRange(-10,10)
        self.p44.setDecimals(4)
        btn4= QtWidgets.QPushButton("LER TENSÃO")

        self.label5 = QtWidgets.QLabel("Coeficiente a", self)
        self.label6 = QtWidgets.QLabel("Coeficiente b", self)
        self.a = QtWidgets.QDoubleSpinBox()
        self.a.setRange(-10000,10000)
        self.a.setDecimals(4)
        self.b = QtWidgets.QDoubleSpinBox()
        self.b.setRange(-10000,10000)
        self.b.setDecimals(4)
        self.button1 = QtWidgets.QPushButton("CALCULAR")
        self.button2 = QtWidgets.QPushButton("AVANÇAR")
        #button3 = QtWidgets.QPushButton("Avançar")
        grid = QtWidgets.QGridLayout()
        grid.setSpacing(4)

        grid.addWidget(self.label1, 1, 0)
        grid.addWidget(alt1, 1, 1)
        grid.addWidget(self.p1, 1, 2)
        grid.addWidget(volt1, 1, 3)
        grid.addWidget(self.p11, 1, 4)
        grid.addWidget(btn1, 1, 5)

        grid.addWidget(self.label2, 2, 0)
        grid.addWidget(alt2, 2, 1)
        grid.addWidget(self.p2, 2, 2)
        grid.addWidget(volt2, 2, 3)
        grid.addWidget(self.p22, 2, 4)
        grid.addWidget(btn2, 2, 5)

        grid.addWidget(self.label3, 3, 0)
        grid.addWidget(alt3, 3, 1)
        grid.addWidget(self.p3, 3, 2)
        grid.addWidget(volt3, 3, 3)
        grid.addWidget(self.p33, 3, 4)
        grid.addWidget(btn3, 3, 5)

        grid.addWidget(self.label4, 4, 0)
        grid.addWidget(alt4, 4, 1)
        grid.addWidget(self.p4, 4, 2)
        grid.addWidget(volt4, 4, 3)
        grid.addWidget(self.p44, 4, 4)
        grid.addWidget(btn4, 4, 5)

        grid.addWidget(self.label5, 5, 0)
        grid.addWidget(self.a, 5, 2,)

        grid.addWidget(self.label6, 6, 0)
        grid.addWidget(self.b, 6, 2,)

        grid.addWidget(self.button1, 7, 0)
        grid.addWidget(self.button2, 7, 1)

        self.instr = QtWidgets.QLabel("Instruções para calibração de sensor: \nPara fazer a calibração do sensor de nível, ligar a bomba no modo manual utilizando o interruptor na bancada, e controlar o\nnível de líquido até os pontos escolhidos.\n\nFaça a leitura visual da altura de líquido na régua do tanque na bancada e preencha o campo 'Altura' do respectivo ponto, depois clique no botão\n'LER TENSÃO' para ler a tensão do sensor no ponto escolhido. Após relizar o procedimento para 4 pontos, clicar em 'CALCULAR' para\ncalcular os coeficentes da equação de calibração do sensor. Após esse procedimento, prossiga clicando em 'AVANÇAR'", self)
        self.instr.setFont(QtGui.QFont('Arial', 10))
        self.instr2 = QtWidgets.QLabel("\n \n \n \n \n \n \nOBS: Recomenda-se a calibração utilizando as alturas mínima e máxima de líquido no tanque, e dois pontos intermediários\nespaçados de forma uniforme", self)
        self.instr2.setFont(QtGui.QFont('Arial', 10))
        grid.addWidget(self.instr, 8, 0, 8, 5)
        grid.addWidget(self.instr2, 9, 0, 9, 5)
        self.setLayout(grid)
        
        
        btn1.clicked.connect(self.ler_volt1)
        btn2.clicked.connect(self.ler_volt2)
        btn3.clicked.connect(self.ler_volt3)
        btn4.clicked.connect(self.ler_volt4)


        self.button1.clicked.connect(self.calc_coef)
        self.button2.clicked.connect(self.avancar)
        

    def avancar(self): # go to main GUI
        CalibrationWindow.a = self.a.value()
        CalibrationWindow.b = self.b.value()
        self.switch_window.emit()


    def a_changed(self):
        CalibrationWindow.a = self.a.value()
    
    def b_changed(self):
        CalibrationWindow.b = self.b.value()

    def calc_coef(self): # calibration calculation
        x1=self.p1.value()
        x2=self.p2.value()
        x3=self.p3.value()
        x4=self.p4.value()
        x=[x1,x2,x3,x4]
        
        y1=self.p11.value()
        y2=self.p22.value()
        y3=self.p33.value()
        y4=self.p44.value()
        y=[y1,y2,y3,y4]

        coefficients = np.polyfit(y,x,1)

        self.a.setValue(coefficients[0])
        self.b.setValue(coefficients[1])

    def ler_volt1(self):
        with nidaqmx.Task() as task:
                task.ai_channels.add_ai_voltage_chan(WorkerThread1.dev+"/"+WorkerThread1.input_h, terminal_config = TerminalConfiguration.RSE) 
                data = task.read(number_of_samples_per_channel=1)
                #print(data[0])
                self.p11.setValue(data[0])
            
    def ler_volt2(self):
        with nidaqmx.Task() as task:
                task.ai_channels.add_ai_voltage_chan(WorkerThread1.dev+"/"+WorkerThread1.input_h, terminal_config = TerminalConfiguration.RSE) 
                data = task.read(number_of_samples_per_channel=1)
                #print(data[0])
                self.p22.setValue(data[0])
    
    def ler_volt3(self):
        with nidaqmx.Task() as task:
                task.ai_channels.add_ai_voltage_chan(WorkerThread1.dev+"/"+WorkerThread1.input_h, terminal_config = TerminalConfiguration.RSE) 
                data = task.read(number_of_samples_per_channel=1)
                #print(data[0])
                self.p33.setValue(data[0])
        
    def ler_volt4(self):
        with nidaqmx.Task() as task:
                task.ai_channels.add_ai_voltage_chan(WorkerThread1.dev+"/"+WorkerThread1.input_h, terminal_config = TerminalConfiguration.RSE) 
                data = task.read(number_of_samples_per_channel=1)
                #print(data[0])
                self.p44.setValue(data[0])


#####################################################################################################

# Controle de troca de interfaces // GUI exchange
class Controller:
    def __init__(self):
        pass

    def show_calibration(self):
        self.calibration = CalibrationWindow()
        self.calibration.switch_window.connect(self.show_main)
        self.calibration.show()

    def show_main(self):
        self.window = MainWindow()
        self.calibration.close()
        self.window.show()
        self.window.open_window.connect(self.show_tuning)

    def show_tuning(self):
        self.tuning = TuningWindow()
        self.tuning.show()

def main():
    app = QtWidgets.QApplication(sys.argv)
    controller = Controller()
    controller.show_calibration()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
