import numpy as np
from numpy import sin, log


from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from scipy.misc import derivative

from tkinter import *
from sympy import *
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

n = -1
dRx = []
dRy = []
dfx = []
x_cur = [[]]
num_points = 50
x_res = []
t1 = []
a = 0
b = 0
def define_diff():
    xa = []
    xb = []
    x = []
    for i in range(n):
        xa.append(Symbol("xa"+str(i + 1)))
        xb.append(Symbol("xb" + str(i + 1)))
        x.append(Symbol("x" + str(i + 1)))
    for i in range(n):
        tmp1 = []
        tmp2 = []
        tmp3 = []
        for k in range(n):
            tmp1.append(diff(edge[i], xa[k]))
            tmp2.append(diff(edge[i], xb[k]))
            tmp3.append(diff(func[i],x[k]))
        dRx.append(tmp1)
        dRy.append(tmp2)
        dfx.append(tmp3)


def solve_ode(x1, x2):
    plot = Tk()
    fig = Figure(figsize=(5, 5), dpi=100)
    plot1 = fig.add_subplot(111)
    plot1.plot(x1, x2)
    plot1.grid()
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
    x_a_p = x_cur[0]
    x_b_p = x_cur[-1]
    dRdx = np.zeros((n, n))
    dRdy = np.zeros((n, n)) #laksjdflkajsefjsaeiojfioasjfiojasef
    xa = []
    xb = []
    for i in range(n):
        xa.append(Symbol("xa"+str(i + 1)))
        xb.append(Symbol("xb" + str(i + 1)))
    for i in range(n):
        for k in range(n):
            tmp_fun = dRx[k][i]
            for j in range(n):
                tmp_fun = tmp_fun.subs(xa[j], x_a_p[j])
                tmp_fun = tmp_fun.subs(xb[j], x_b_p[j])
            dRdx[k][i] = tmp_fun.evalf()
    for i in range(n):
        for k in range(n):
            tmp_fun = dRy[k][i]
            for j in range(n):
                tmp_fun = tmp_fun.subs(xa[j], x_a_p[j])
                tmp_fun = tmp_fun.subs(xb[j], x_b_p[j])
            dRdy[k][i] = tmp_fun.evalf()
    return [dRdx, dRdy]

def Xrhs(tt, y):
    t = Symbol("t")
    x = []
    x_c = x_cur[ind(tt, a, b)]
    for i in range(1, 7):
        x.append(Symbol("x" + str(i)))
    dfdxp = np.zeros((n, n))
    for i in range(n):
        for k in range(n):
            tmp_fun = dfx[k][i]
            tmp_fun.subs(t, tt)
            for j in range(n):
                tmp_fun = tmp_fun.subs(x[j], x_c[j])
            dfdxp[k][i] = tmp_fun.evalf()
    return np.matmul(dfdxp, np.array(y))

def X(p, a, b):
    res = []
    t = np.linspace(a, b, num=5)
    for i in range(n):
        init = np.zeros(n)
        init[i] = 1
        res.append(solve_ivp(Xrhs, [a, b], init).y[:,-1])
    return np.array(res).transpose()


def ind(t, a, b):
    res = int(np.floor((t - a) * (num_points - 1) / (b - a)))
    if (res < 0):
        res = 0
    if (res > num_points - 1):
        res = num_points - 1
    return res


def dfidmu(p, a, b, t0):
    t = np.linspace(a, b, num=num_points)
    global x_cur
    x_cur = odeint(fun_mat, p, np.flip(t[0:(ind(t0, a, b) + 1)]))
    x_cur = np.flip(x_cur, 0)
    x_cur = np.append(x_cur, odeint(fun_mat, p, t[ind(t0, a, b):])[1:], axis=0)

    R = dRdxdy(p, a, b, t0)
    Xa = X(p, t0, a)
    Xb = X(p, t0, b)
    res = np.matmul(R[0], Xa) + np.matmul(R[1], Xb)
    return res

def prhs(mu, p, a, b, t0, f0):

    print(p, " ", mu)
    pb['value'] = mu * 100
    window.update()
    f_inv = np.linalg.inv(dfidmu(p, a, b, t0))
    return -np.matmul(f_inv, f0)

def test():
    func.append(parse_expr("x3*x2"))
    func.append(parse_expr("x3*(-x1+sin(x2))"))
    func.append(parse_expr("0"))
    func.append(parse_expr("0"))
    edge.append(parse_expr("xa1-xa4"))
    edge.append(parse_expr("xb1-xb4"))
    edge.append(parse_expr("xa2"))
    edge.append(parse_expr("xb2"))
    global a, b, n
    a = 0
    b = 1
    n = 4


pb = 0
wait = 0

def draw():
    global t1, x_res
    conv = {
        "t" : t1
    }
    for i in range(n):
        conv["x" + str(i + 1)] = x_res[:, i]
    solve_ode(conv[dr1.get()], conv[dr2.get()])


def tm():
    window.update()
    global a, b
    global x_res
    global t1
    global func, edge
    func = []
    edge = []
    button_draw["state"] = DISABLED
    window.update()
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


    #test()
    #p0 = [2,0,2*3.1415926,2]
    #t0 = 0

    define_diff()
    f0 = fi0(p0, a, b, t0)
    rhs = lambda tt, y: prhs(tt, y, a, b, t0, f0)
    mu = np.linspace(0, 1, num=3)
    p_res = solve_ivp(rhs, [0,1], p0, method='RK23').y[:,-1]
    t1 = np.linspace(t0, a)
    x_res = odeint(fun_mat, p_res, t1)
    x_res = np.flip(x_res, 0)
    t = np.linspace(t0, b)
    x_res = np.append(x_res, odeint(fun_mat, p_res, t), axis=0)
    t1 = np.flip(t1)
    t1 = np.append(t1, t)
    button_draw["state"] = NORMAL
    window.update()



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

fr1 = Frame(fr)
fr1.pack(side=TOP)
Label(fr1, text="Оси графиков:", font=font).pack(side="left", expand=True, pady=5, padx=15)

dr1 = ttk.Combobox(
    fr1,
    width = 5,
    state="readonly",
    values=["t", "x1", "x2", "x3", "x4", "x5", "x6"]
)
dr2 = ttk.Combobox(
    fr1,
    width = 5,
    state="readonly",
    values=["t", "x1", "x2", "x3", "x4", "x5", "x6"]
)
dr1.pack(pady=16, side=LEFT)
dr2.pack(pady=16, side=LEFT)
button_draw = Button(fr1, text="Нарисовать", command=draw, state=DISABLED)
button_draw.pack(side=LEFT, fill=BOTH, pady=16, padx=15)


fr = Frame(window)
fr.pack(side=BOTTOM)
Button(fr, text="Выполнить", command=tm).pack(side=BOTTOM, fill=BOTH, pady=16)
pb = ttk.Progressbar(window, orient="horizontal", length=200, value=0)
pb.pack(pady=16, side=BOTTOM)


func = []
edge = []

window.mainloop()


