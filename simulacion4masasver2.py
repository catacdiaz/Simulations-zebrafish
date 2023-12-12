import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

#parámetros

k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = 10, 10, 10, 10, 10, 10, 1, 1, 1, 1
L1, L2, L3, L4, L5, L6, L7, L8, L9, L10 = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1

b = 3   #coeficiente de amortiguamiento
x0 = 1   #espacio entre resortes
F2 = 0.2 #fuerza externa
F3 = 0 #fuerza externa
F4 = 0.4 #fuerza externa

# Ecuaciones de movimiento
def ecuaciones_movimiento(t, y):
    x1, y1, x2, y2, x3, y3, x4, y4 = y

    # Longitudes y fuerzas del resorte
    r1 = np.sqrt(x1**2 + y1**2)
    r2 = np.sqrt((x2-x0)**2 + y2**2)
    r3 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    r4 = np.sqrt((x3-x1)**2 + (y3-y1)**2)
    r5 = np.sqrt((x4-x2)**2 + (y4-y2)**2)
    r6 = np.sqrt((x4-x3)**2 + (y4-y3)**2)
    r7 = np.sqrt((x2)**2 + (y2)**2)
    r8 = np.sqrt((x0-x1)**2 + (y1)**2)
    r9 = np.sqrt((x4-x1)**2 + (y4-y1)**2)
    r10 = np.sqrt((x2-x3)**2 + (y3-y2)**2)

    dx1dt  = -k1 * x1 * (r1 - L1) /(r1*b) + k3 * (x2 - x1)*(r3 - L3)/(r3*b) + k4*(x3-x1)*(r4 - L4)/(r4*b)  +k8*(x0-x1)*(r8 - L8)/(r8*b) + k9*(x4-x1)*(r9 - L9)/(r9*b)
    dy1dt = -k1 * y1 * (r1 - L1) /(r1*b) + k3 * (y2 - y1)*(r3 - L3)/(r3*b) + k4*(y3-y1)*(r4 - L4)/(r4*b)   -k8*(y1)*(r8 - L8)/(r8*b) + k9*(y4-y1)*(r9 - L9)/(r9*b)

    dx2dt  = -k2 * (x2 - x0) * (r2 - L2)/(r2*b) - k3*(x2 - x1)*(r3 - L3)/(r3*b) + k5*(x4-x2)*(r5 - L4)/(r5*b)   -k7*(x2)*(r7 - L7)/(r7*b) - k10*(x2-x3)*(r10 - L10)/(r10*b) + F2
    dy2dt = -k2 * y2 * (r2 - L2) /(r2*b) - k3 * (y2 - y1)*(r3 - L3)/(r3*b) + k5*(y4-y2)*(r5 - L4)/(r5*b)   -k7*(y2)*(r7 - L7)/(r7*b) + k10*(y2-y3)*(r10 - L10)/(r10*b)

    dx3dt  = -k4 * (x3-x1) * (r4 - L4) /(r4*b) + k6*(x4 - x3)*(r6 - L6)/(r6*b)   + k10*(x2-x3)*(r10 - L10)/(r10*b) + F3
    dy3dt = -k4 * (y3-y1) * (r4 - L4) /(r4*b) + k6*(y4 - y3)*(r6 - L6)/(r6*b)   -  k10*(y2-y3)*(r10 - L10)/(r10*b)

    dx4dt  = -k5 * (x4-x2) * (r4 - L4) /(r4*b) - k6*(x4 - x3)*(r6 - L6)/(r6*b)  - k9*(x4-x1)*(r9 - L9)/(r9*b) + F4
    dy4dt = -k5 * (y4-y2) * (r4 - L4) /(r4*b) - k6*(y4 - y3)*(r6 - L6)/(r6*b)  - k9*(y4-y1)*(r9 - L9)/(r9*b)


    return [dx1dt, dy1dt, dx2dt, dy2dt, dx3dt, dy3dt, dx4dt, dy4dt]

#condiciones iniciales
x1_0 = 0 
y1_0 = 1 

x2_0 = 1 
y2_0 = 1  

x3_0 = 0 
y3_0 = 2  

x4_0 = 1
y4_0 = 2 

#tiempo
t_span = (0, 3)
t_eval = np.linspace(t_span[0], t_span[1], 100)

#condiciones iniciales
y0 = [x1_0, y1_0, x2_0, y2_0, x3_0, y3_0, x4_0, y4_0]


sol = solve_ivp(ecuaciones_movimiento, t_span, y0, t_eval=t_eval, method='RK45')

#renombre soluciones
x1 = sol.y[0]
y1 = sol.y[1]
x2 = sol.y[2]
y2 = sol.y[3]
x3 = sol.y[4]
y3 = sol.y[5]
x4 = sol.y[6]
y4 = sol.y[7]

#cánculo de ángulos
pi = np.array([np.pi for _ in range(len(sol.t))])
x_0 = np.array([x0 for _ in range(len(sol.t))])

tan1 = x1/y1
tan2 = (x2 - x_0)/y2
tan3 = x3/y3
tan4 = (x4 - x_0)/y4

theta1 = np.arctan(tan1)*180/np.pi
theta2 = np.arctan(tan2)*180/np.pi
theta3 = np.arctan(tan3)*180/np.pi
theta4 = np.arctan(tan4)*180/np.pi

theta_3 = (theta3 + 90)
theta_4 = (90 - theta4) 

f = (np.abs(theta2) - np.abs(theta1))/(np.abs(theta2) + np.abs(theta1))


alpha = (90 - theta2)
betha = (90 + theta1)

f_1 = (-theta2 -theta1)/(-theta2 +theta1 + 180)
f_2 = (-theta4 -theta3)/(-theta4 +theta3 + 180)


# Animación
fig, ax = plt.subplots()
dot1 = ax.scatter([], [], s=20, label='Masa 1') 
dot2 = ax.scatter([], [], s=20, label='Masa 2')
dot3 = ax.scatter([], [], s=20, label='Masa 3')  
dot4 = ax.scatter([], [], s=20, label='Masa 4')


#lineas
line1, = ax.plot([0, 0], [0, 0], 'r--', alpha=0.5)
line2, = ax.plot([1, 1], [0, 0], 'b--', alpha=0.5)
line3, = ax.plot([0, 0], [0, 0], 'g--', alpha=0.5)
line4, = ax.plot([0, 0], [0, 0], 'p--', alpha=0.5)
line5, = ax.plot([0, 0], [0, 0], 'g--', alpha=0.5)
line6, = ax.plot([0, 0], [0, 0], 'p--', alpha=0.5)
line7, = ax.plot([0, 0], [0, 0], 'r--', alpha=0.5)
line8, = ax.plot([1, 1], [0, 0], 'b--', alpha=0.5)
line9, = ax.plot([0, 0], [0, 0], 'g--', alpha=0.5)
line10, = ax.plot([0, 0], [0, 0], 'p--', alpha=0.5)

ax.set_xlim(0, 1.5)
ax.set_ylim(0, 2.5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Movimiento de cuatro masas con resortes')
ax.legend()

def init():
    dot1.set_offsets(np.empty((0, 2)))
    dot2.set_offsets(np.empty((0, 2)))
    dot3.set_offsets(np.empty((0, 2)))
    dot4.set_offsets(np.empty((0, 2)))
    return dot1, dot2, dot3, dot4, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10

def update(frame):
    x1 = sol.y[0][frame]
    y1 = sol.y[1][frame]
    x2 = sol.y[2][frame]
    y2 = sol.y[3][frame]
    x3 = sol.y[4][frame]
    y3 = sol.y[5][frame]
    x4 = sol.y[6][frame]
    y4 = sol.y[7][frame]

    dot1.set_offsets([[x1, y1]])
    dot2.set_offsets([[x2, y2]])
    dot3.set_offsets([[x3, y3]])
    dot4.set_offsets([[x4, y4]])

    line1.set_data([0, x1], [0, y1])
    line2.set_data([1, x2], [0, y2])
    line3.set_data([x1, x2], [y1, y2])
    line4.set_data([x1, x3], [y1, y3])
    line5.set_data([x2, x4], [y2, y4])
    line6.set_data([x3, x4], [y3, y4])
    line7.set_data([0, x2], [0, y2])
    line8.set_data([1, x1], [0, y1])
    line9.set_data([x1, x4], [y1, y4])
    line10.set_data([x2, x3], [y2, y3])     

    return dot1, dot2, dot3, dot4, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10

ani = FuncAnimation(fig, update, frames=len(sol.t), init_func=init, blit=False)

plt.show()


plt.plot(sol.t, theta1, label='theta1')
plt.plot(sol.t, theta2, label='theta2')
plt.xlabel('Tiempo')
plt.ylabel('Ángulo')
plt.title('Ángulos')
plt.legend()
plt.show()


plt.plot(sol.t, betha, label='b')
plt.plot(sol.t, alpha, label='a')
plt.xlabel('Tiempo')
plt.ylabel('Ángulo')
plt.title('Ángulos con respecto a horizontal')
plt.legend()
plt.show()

plt.plot(sol.t, betha, label='b')
plt.plot(sol.t, alpha, label='a')
plt.xlabel('Tiempo')
plt.ylabel('Ángulo')
plt.title('Ángulos con respecto a horizontal')
plt.legend()
plt.show()

plt.plot(sol.t, theta_3, label='theta3')
plt.plot(sol.t, theta_4, label='theta4')
plt.xlabel('Tiempo')
plt.ylabel('Ángulo')
plt.title('Ángulos con respecto a horizontal')
plt.legend()
plt.show()

plt.plot(f_1, sol.t, label='f')
plt.xlabel('Ángulo')
plt.ylabel('Tiempo')
plt.title(f'Relación f con ángulos alpha y betha,con m1 y m2')
plt.legend()
plt.show()

plt.plot(f_2, sol.t, label='f')
plt.xlabel('Ángulo')
plt.ylabel('Tiempo')
plt.title(f'Relación f con ángulos alpha y betha,con m3 y m4')
plt.legend()
plt.show()

"""
plt.plot(sol.t, sol.y[0], label='x1')
plt.plot(sol.t, sol.y[1], label='y1')
plt.xlabel('Tiempo')
plt.ylabel('Posición')
plt.title('Masa 1')
plt.legend()
plt.show()


plt.plot(sol.t, sol.y[2], label='x2')
plt.plot(sol.t, sol.y[3], label='y2')
plt.xlabel('Tiempo')
plt.ylabel('Posición')
plt.title('Masa 2')
plt.legend()
plt.show()



plt.plot(sol.t, sol.y[4], label='x3')
plt.plot(sol.t, sol.y[5], label='y3')
plt.xlabel('Tiempo')
plt.ylabel('Posición')
plt.title('Masa 3')
plt.legend()
plt.show()


plt.plot(sol.t, sol.y[6], label='x4')
plt.plot(sol.t, sol.y[7], label='y4')
plt.xlabel('Tiempo')
plt.ylabel('Posición')
plt.title('Masa 4')
plt.legend()
plt.show()
"""