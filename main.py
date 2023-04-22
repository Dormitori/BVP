import numpy as np
from numpy import sin, log



from scipy.integrate import odeint
from scipy.misc import derivative

from tkinter import *
from sympy import *

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

n = -1

def solve_ode(y, t):
    plot = Tk()
    fig = Figure(figsize=(5, 5), dpi=100)
    plot1 = fig.add_subplot(111)
    plot1.plot(t, y)

    canvas = FigureCanvasTkAgg(fig, master=plot)
    canvas.draw()

    canvas.get_tk_widget().pack()

    toolbar = NavigationToolbar2Tk(canvas, plot)
    toolbar.update()

    # placing the toolbar on the Tkinter window
    canvas.get_tk_widget().pack()

    plot.mainloop()


def fun_mat(y, t):
    tt = Symbol("t")
    x = []
    for i in range(1, 7):
        x.append(Symbol("x" + str(i)))
    global n
    res = np.zeros(n)
    for i in range(n):
        tmp_func = func[i].subs(tt, t)
        for j in range(n):
            tmp_func = tmp_func.subs(x[j],y[j])
        res[i] = tmp_func.evalf()
    return res

def fi0(p0, a, b, t0):
    t = np.linspace(t0, a)
    res_ode = odeint(fun_mat, p0, t)
    x_a_p0 = res_ode[-1]
    t = np.linspace(t0, b)
    res_ode = odeint(fun_mat, p0, t)
    x_b_p0 = res_ode[-1]

    xa = []
    xb = []
    for i in range(1, 7):
        xa.append(Symbol("xa" + str(i)))
        xb.append(Symbol("xb" + str(i)))
    res = np.zeros(n)
    for i in range(n):
        tmp_func = edge[i]
        for j in range(n):
            tmp_func = tmp_func.subs(xa[j], x_a_p0[j])
            tmp_func = tmp_func.subs(xb[j], x_b_p0[j])
        res[i] = tmp_func.evalf()
    return res

def dRdxdy(p, a, b, t0):
    t = np.linspace(t0, a)
    res_ode = odeint(fun_mat, p, t)
    x_a_p = res_ode[-1]
    t = np.linspace(t0, b)
    res_ode = odeint(fun_mat, p, t)
    x_b_p = res_ode[-1]
    dRdx = np.zeros((n, n))
    dRdy = np.zeros((n, n)) #laksjdflkajsefjsaeiojfioasjfiojasef
    xa = []
    xb = []
    for i in range(n):
        xa.append(Symbol("xa"+str(i + 1)))
        xb.append(Symbol("xb" + str(i + 1)))
    for i in range(n):
        for k in range(n):
            tmp_fun = diff(edge[k], xa[i])
            for j in range(n):
                tmp_fun = tmp_fun.subs(xa[j], x_a_p[j])
                tmp_fun = tmp_fun.subs(xb[j], x_b_p[j])
            dRdx[k][i] = tmp_fun.evalf()
    for i in range(n):
        for k in range(n):
            tmp_fun = diff(edge[k], xb[i])
            for j in range(n):
                tmp_fun = tmp_fun.subs(xa[j], x_a_p[j])
                tmp_fun = tmp_fun.subs(xb[j], x_b_p[j])
            dRdy[k][i] = tmp_fun.evalf()
    return [dRdx, dRdy]

def Xrhs(y, tt, p):
    t = Symbol("t")
    x = []
    for i in range(1, 7):
        x.append(Symbol("x" + str(i)))
    dfdxp = np.zeros((n, n))
    for i in range(n):
        for k in range(n):
            tmp_fun = diff(func[k], x[i])
            tmp_fun.subs(t, tt)
            for j in range(n):
                tmp_fun = tmp_fun.subs(x[j], p[j])
            dfdxp[k][i] = tmp_fun.evalf()
    return np.matmul(dfdxp, np.array(y))

def X(p, a, b):
    rhs = lambda y, tt: Xrhs(y, tt, p)
    res = []
    t = np.linspace(a, b)
    for i in range(n):
        init = np.zeros(n)
        init[i] = 1
        res.append(odeint(rhs, init, t)[-1])
    return np.array(res).transpose()

def dfidmu(p, a, b, t0):
    R = dRdxdy(p, a, b, t0)
    Xa = X(p, t0, a)
    Xb = X(p, t0, b)
    res = np.matmul(R[0], Xa) + np.matmul(R[1], Xb)
    return res

def prhs(p, mu, a, b, t0, f0):
    f_inv = np.linalg.inv(dfidmu(p, a, b, t0))
    return -np.matmul(f_inv, f0)

def tm():
    global func, edge
    func = []
    edge = []
    for i in range(1, 7):
        input = texts_diff[i - 1].get("1.0", "end-1c")
        if input != "":
            func.append(parse_expr(input, evaluate=True))
    for i in range(6):
        input = texts_edge[i].get("1.0", "end-1c")
        if input != "":
            edge.append(parse_expr(input, evaluate=True))
    global n
    print(n)
    n = len(func)
    print(n)
    a = eval(init_a.get("1.0", "end-1c"))
    b = eval(init_b.get("1.0", "end-1c"))
    p0 = np.array(eval(init_value.get("1.0", "end-1c")))
    t0 = eval(init_time.get("1.0", "end-1c"))
    f0 = fi0(p0, a, b, t0)
    rhs = lambda y, tt: prhs(y, tt, a, b, t0, f0)
    mu = np.linspace(0, 1)
    p_res = odeint(rhs, p0, mu)[-1]
    t1 = np.linspace(t0, a)
    x_res = odeint(fun_mat, p_res, t1)
    x_res = np.flip(x_res)
    t = np.linspace(t0, b)
    x_res = np.append(x_res, odeint(fun_mat, p_res, t), axis=0)
    t1 = np.flip(t1)
    t1 = np.append(t1, t)
    print(x_res)
    solve_ode(x_res, t1)


window = Tk()
window.title("BVP")
window.geometry("1200x500")


font = ("Comic Sans MS", 15)
texts_diff = []


fr = Frame(window)
fr.pack(side=LEFT)

fr1 = Frame(fr)
fr1.pack(side=TOP)
Label(fr1, text="Введите систему дифф уравнений:", font=font).pack(side="top", expand=True, pady=5)

for i in range(1, 7):
    fr1 = Frame(fr)
    fr1.pack(side=TOP)
    texts_diff.append(Text(fr1, width=30, height=1, font=("Arial", 15)))
    Label(fr1, text="dx" + str(i) + "/dt: ", font=font).pack(side="left", expand=True, pady=5, padx=40)
    texts_diff[i - 1].pack(side="left", fill="x", expand=True)

fr1 = Frame(fr)
fr1.pack(side=TOP)
Label(fr1, text="Временной промежуток:", font=font).pack(side="left", expand=True, pady=5)
init_a = Text(fr1, width=5, height=1, font=("Arial", 15))
init_a.pack(side="left", expand=True, pady=5, padx=10)
init_b = Text(fr1, width=5, height=1, font=("Arial", 15))
init_b.pack(side="left", expand=True, pady=5)

fr1 = Frame(fr)
fr1.pack(side=TOP)
Label(fr1, text="Начальное время:", font=font).pack(side="left", expand=True, pady=5)
init_time = Text(fr1, width=20, height=1, font=("Arial", 15))
init_time.pack(side="left", expand=True, pady=5)


fr = Frame(window)
fr.pack(side=RIGHT)

fr1 = Frame(fr)
fr1.pack(side=TOP)
Label(fr1, text="Введите систему начальных условий:", font=font).pack(side="top", expand=True, pady=5)

texts_edge = []

for i in range(1, 7):
    fr1 = Frame(fr)
    fr1.pack(side=TOP)
    texts_edge.append(Text(fr1, width=30, height=1, font=("Arial", 15)))
    texts_edge[i - 1].pack(side="left", fill="x", expand=True)
    Label(fr1, text=" = 0", font=font).pack(side="left", expand=True, pady=5, padx=40)



fr1 = Frame(fr)
fr1.pack(side=TOP)
Label(fr1, text="Вектор начального приближения:", font=font).pack(side="left", expand=True, pady=5)
init_value = Text(fr1, width=10, height=1, font=("Arial", 15))
init_value.pack(side="left", expand=True, pady=5, padx=40)



fr = Frame(window)
fr.pack(side=BOTTOM)
Button(window, text="Выполнить", command=tm).pack(side=BOTTOM,fill=BOTH, pady=16)

func = []
edge = []

window.mainloop()


