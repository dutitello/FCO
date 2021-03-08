# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:44:35 2019

Programa para verificação de seções de concreto armado de formato genérico
para flexo-compressão oblíqua de acordo com a NBR-6118:2014

@author: MVREAL
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#    
#
def tensao(es1,esm,fyd):
#
#
#  Calcula a tensão no aço
#  es = módulo de elasticidade do aço em kN/cm2
#  es1 = deformação de entrada
#  fyd = tensão de escoamento de cálculo em kN/cm2
#  ss1 = tensão de saída em kN/cm2
#
#  Trabalhando com deformação positiva
    ess = np.abs(es1);
    eyd = fyd / esm;
    if ess < eyd:
#
        sigmasd = esm * ess;
#
    else:
#
        sigmasd = fyd;
#
#  Trocando o sinal se necessário
    if es1 < 0:
#
        sigmasd = -sigmasd;
    return sigmasd
#
#        
def funcao(nsd,qsi,alfa,nc1,xc,yc,alfac,sigmacd,ns,xs,ys,eu,e0,alamb,esm,fyd,astotal):       
#
# Calcula o valor da função f(x0)dada na equação (2.5.11) do Volume 3 de Curso de Concreto Armado
#
#  qsi=x0/h é a profundidade relativa da linha neutra
#  rc é a resultante de compressão do concreto adimensional dada na equação (2.4.4)
#  bc é a posição da resultante adimensional dada na equação (2.4.5)
#  soma1 é o somatório contido no denominador da equação (2.5.9)
#  soma2 é o somatatório contido no denominador da equação (2.5.10)
#  f é o resultado da equação (2.5.11)
#
# Todas as variáveis de entrada são públicas     
# Os parâmetros rc,bc,soma1,soma2,f são calculados nessa sub-rotina
#
#
# Cálculo das coordenadas dos vértices da poligonal segundo o sistema x'oy'
# paralelo e normal a linha neutra com orientação alfa em relação ao eixo x.
#
    xl=np.zeros(nc1)
    yl=np.zeros(nc1)
#   
    for i in range(nc1):
        xl[i]=xc[i]*np.cos(alfa)+yc[i]*np.sin(alfa)
        yl[i]=yc[i]*np.cos(alfa)-xc[i]*np.sin(alfa)
#
#    print(xl,yl)
# Valores máximo e mínimo de y'
# Cálculo da altura h da seção
#
    ylmax=np.amax(yl)
    ylmin=np.amin(yl)
    h=np.abs(ylmax-ylmin)
#    print("h =",h)
#

# Cálculo das coordenadas das armaduras segundo o sistema x'oy'
# paralelo e normal a linha neutra com orientação alfa em relação ao eixo x.
# Montagem do vetor Beta(i)=dsi/h
#
    xsl=np.zeros(ns)
    ysl=np.zeros(ns)
    beta=np.zeros(ns)
    for i in range(ns):
        xsl[i]=xs[i]*np.cos(alfa)+ys[i]*np.sin(alfa)
        ysl[i]=ys[i]*np.cos(alfa)-xs[i]*np.sin(alfa)
        beta[i]=np.abs(ylmax-ysl[i])/h
#
#    print(xsl,ysl)
#    print(beta)
# Valores máximo e mínimo de ys'
# Cálculo da altura útil d  da seção
#
    yslmin=np.amin(ysl)
    d=np.abs(ylmax-yslmin)
#    print("d =",d)

    betamax=d/h
#
# Cálculo dos parâmetros qsi e delta
#
    x0=qsi*h
#    print("x0 =",x0)
#           
#  Constantes para o cálculo das deformações das camadas da armadura
#  Observar que o primeiro índice é zero
    q1 = eu * betamax / (eu + 10);
    if qsi <= q1:
#
#  A linha neutra está no domínio 2
        c = 0.01 / (betamax - qsi);
#
    elif qsi<=1:
#            
#  A linha neutra está nos domínios 3,4 e 4a
        c = eu / (1000 * qsi);
#
    else:
#        
#  A linha neutra está no domínio 5
        c = (e0 / 1000) / (qsi - akapa);
       
#  
#    print("c =",c)  
#        
# Resultante de compressão no concreto
#
#
    if alamb*x0<h:
#
# Seção parcialmente comprimida
#
#    
# Coordenada ycl do início da zona de concreto comprimido
# com o diagrama retangular simplificado
#
        ycl = ylmax-alamb*x0
#        print("ycl =",ycl)
        kint=0
        xaux=np.zeros(100)
        yaux=np.zeros(100)
        xaux[0:nc1]=xl[0:nc1]
        yaux[0:nc1]=yl[0:nc1]
        for j in range(nc1-1):
            deltaxl=xl[j+1]-xl[j]
            deltayl=yl[j+1]-yl[j]
#            print("j =",j,"deltaxl =",deltaxl,"deltayl =",deltayl)
            if np.abs(deltayl)>=10e-8:
                xint=xl[j]+(ycl-yl[j])*deltaxl/deltayl
                if((xint>=xl[j])and(xint<=xl[j+1])or((xint<=xl[j])and(xint>=xl[j+1]))):
#                    print("xint =",xint,"j =",j,xl[j],xl[j+1])
                    kint+=1
#                    print("kint =",kint)
                    nv=nc1+kint
                    xaux[j+kint]=xint
                    yaux[j+kint]=ycl
                    xaux[j+kint+1:nv]=xl[j+1:nc1]
                    yaux[j+kint+1:nv]=yl[j+1:nc1]
                nv=nc1+kint
#        print("nv =",nv)
#        print("aux",xaux[0:nv],yaux[0:nv])
        xl=np.zeros(nv)
        yl=np.zeros(nv)
        xl[0:nv]=xaux[0:nv]
        yl[0:nv]=yaux[0:nv]
#
#        print("xl =",xl,"yl =",yl)
# Determinação dos vértices da zona de concreto comprimido
#
        kvc=0
        xccl=np.zeros(nv)
        yccl=np.zeros(nv)
        for i in range(nv):
            if yl[i]>=ycl:
                xccl[kvc]=xl[i]
                yccl[kvc]=yl[i]
                kvc=kvc+1
        xccl[kvc]=xccl[0]
        yccl[kvc]=yccl[0]
#
# Seção totalmente comprimida
#
#
    else:
        nv=nc1
        kvc=nc1-1
        nvc=kvc
        xccl=np.zeros(nv)
        yccl=np.zeros(nv)
        xccl[0:nv]=xl[0:nv]
        yccl[0:nv]=yl[0:nv]
#
#    print("xccl =",xccl,"yccl =",yccl)
    nvc=kvc
    xcc=np.zeros(nvc+1)
    ycc=np.zeros(nvc+1)
    for i in range(nvc+1):
        xcc[i]=xccl[i]*np.cos(alfa)-yccl[i]*np.sin(alfa)
        ycc[i]=yccl[i]*np.cos(alfa)+xccl[i]*np.sin(alfa)
    acc=0.00
    sxc=0.00
    syc=0.00
    xcc[nvc]=xcc[0]
    ycc[nvc]=ycc[0]
#
#    print("nv =",nvc)
#    print("xcc =",xcc,"ycc =",ycc)
#
    for k in range(nvc):
        deltaxcc=xcc[k+1]-xcc[k]
        deltaycc=ycc[k+1]-ycc[k]
        acc=acc+1./2.*deltaycc*(2.*xcc[k]+deltaxcc)
        sxc=sxc+1./6.*deltaycc*(3.*xcc[k]*xcc[k+1]+deltaxcc**2)
        syc=syc-1./6.*deltaxcc*(3.*ycc[k]*ycc[k+1]+deltaycc**2)
#
#    print("Acc =",acc)
#    print("Sxc =",sxc)
#    print("Syc =",syc)
#
# Cálculo das tensões nas armaduras
#
    somas = 0;
    somasx = 0;
    somasy = 0;
#
    sigmas=np.zeros(ns)    
    for i in range(ns):
#
        es1 = c * (qsi - beta[i]);
        sigmas[i]=tensao(es1,esm,fyd);
        somas = somas + astotal/ns * sigmas[i];
        somasx = somasx + astotal/ns * xs[i] * sigmas[i];
        somasy = somasy + astotal/ns * ys[i] * sigmas[i];
#
#    print("ds[i]=",beta*h)
#    print("Sigmas =",sigmas)
#    print(somas,somasx,somasy)
    
#
# Esforços resultantes Nrd, Mxrd, Myrd
#
    nrd=sigmacd*acc+somas
    mxrd=sigmacd*sxc+somasx
    myrd=sigmacd*syc+somasy
#
#  Funcao f(x)
#
    f = nsd-sigmacd*acc-somas
#    print(f,nrd,mxrd,myrd)
    return f,nrd,mxrd,myrd
#
#
#
# Entrada de dados
#
# Lê dados de uma planilha do Excel
#
with open('dados.xlsx', 'rb') as target:
    sheet =  pd.read_excel(target, sheet_name='Planilha1')
    data  =  sheet.values
#
nc =np.int( data[0,0] )     # coluna A: nc
nc1=nc+1
xc=np.zeros(nc1)
yc=np.zeros(nc1)
xc[0:nc1] =data[0:nc1,1]    # coluna B: xc
yc[0:nc1] =data[0:nc1,2]     # coluna C: yc
astotal =np.float(data[0,3]) # coluna D: Astotal
ns =np.int(data[0,4])      # coluna E: ns
xs=np.zeros(ns)
ys=np.zeros(ns)
xs =data[0:ns,5]     # coluna F: xs
ys =data[0:ns,6]     # coluna G: ys
fck =np.float(data[0,7])      # coluna H: fck
fyk =np.float(data[0,8])      # coluna I: fyk
esm =np.float(data[0,9])     # coluna J: Es
gamac =np.float( data[0,10])  # coluna K: gamac
gamas =np.float( data[0,11])  # coluna L: gamas
gamaf =np.float( data[0,12]) # coluna M: gamaf
nsk =np.float( data[0,13])      # coluna N: nsk
mxsk =np.float( data[0,14])     # coluna O: ex
mysk =np.float( data[0,15])     # coluna N: ey
alfa0=np.float(data[0,16])      # coluna Q: alfa0
alfaf=np.float(data[0,17])      # coluna R: alfaf
#
#
# Final da entrada de dados
#
#
# Início dos cálculos
#
# Parâmetros do diagrama retangular
if (fck <= 50):
    alamb = 0.8;
    alfac = 0.85;
    eu = 3.5;
    e0 = 2;
else:
    alamb = 0.8 - (fck - 50) / 400;
    alfac = 0.85 * (1 - (fck - 50) / 200);
    a = (90 - fck) / 100;
    eu = 2.6 + 35 * a**4;
    aux = fck - 50;
    a = 0.53;
    e0 = 2 + 0.085 *aux**a;
# Parâmetro kapa que define o ponto com deformação igual a eo no domínio 5
akapa = 1 - e0 / eu;
#            
# Conversão de unidades: transformando para kN e cm
fck = fck / 10;
fyk = fyk / 10;
esm = 100 * esm;
#           
# Resistências de cálculo
fcd = fck / gamac;
sigmacd =0.95*alfac * fcd;
fyd = fyk / gamas;
#
# Esforços solicitantes de cálculo
nsd = gamaf * nsk;
mxsd = gamaf * mxsk;
mysd=gamaf*mysk;
#
# Cálculo das coordenadas do centroide da seção de concreto
# Mudança sistema de referência para o centroide da seção de concreto
#
ac = 0.00
sx = 0.00
sy = 0.00
# 
xc[nc]=xc[0]
yc[nc]=yc[0]
for i in range(nc):  
    xi = xc[i]
    yi = yc[i]
    dx = xc[i + 1] - xc[i]
    dy = yc[i + 1] - yc[i]
#
    ac = ac + (xi + dx / 2.00) * dy
    sx = sx + (xi * (yi + dy / 2.00) + dx * (yi / 2.00 + dy / 3.00)) * dy
    sy = sy + (xi * (xi + dx) + dx ** 2 / 3.00) * dy / 2.00
#
xg = sy / ac
yg = sx / ac
nc1=nc+1
for i in range(nc1):
    xc[i]=xc[i]-xg
    yc[i]=yc[i]-yg
#
for i in range(ns):
    xs[i]=xs[i]-xg
    ys[i]=ys[i]-yg
#
# Coordenadas atualizadas para o centroide
#
#
# Conversão dos ângulos alfa0 e alfaf para radianos
#
nalfa=400
nalfa1=nalfa+1
alfa0=np.pi/180*alfa0
alfaf=np.pi/180*alfaf
deltaalfa=np.abs(alfaf-alfa0)/nalfa
alfa=-deltaalfa
theta=np.zeros(nalfa1)
nid=np.zeros(nalfa1)
mixd=np.zeros(nalfa1)
miyd=np.zeros(nalfa1)
k=-1
#
# Variação da inclinação da linha neutra alfa entre alfa0 e alfaf
#
while (k<nalfa):
    alfa=alfa+deltaalfa
    k=k+1
#    print("k =",k)
#    print("alfa =",alfa*180.00/np.pi)
#
#            
#  Processo iterativo da bissecante
#  
#  Determinação do intervalo solução
#            
#  Valor inicial para a linha neutra adminesional qsi=x/h
    qi = 0;
#  Chamar sub-rotina para calcular o valor da função fi=f(qi)
    qsi = qi;
    f,nrd,mxrd,myrd=funcao(nsd,qsi,alfa,nc1,xc,yc,alfac,sigmacd,ns,xs,ys,eu,e0,alamb,esm,fyd,astotal)
    fi=f
#
#  Valor final para a linha neutra adimensional qsi=x/h
    qf = 1000;
#  Chamar sub-rotina para calcular o valor da função ff=f(qf)
    qsi = qf;
    f,nrd,mxrd,myrd=funcao(nsd,qsi,alfa,nc1,xc,yc,alfac,sigmacd,ns,xs,ys,eu,e0,alamb,esm,fyd,astotal)
    ff=f
    prod = fi * ff;
#  Modificando os extremos do intervalo solução até que prod<=0
    while prod > 0:
#            
        qi = qf;
        fi = ff;
        qf = 10 * qf;
        qsi = qf;
        f,nrd,mxrd,myrd=funcao(nsd,qsi,alfa,nc1,xc,yc,alfac,sigmacd,ns,xs,ys,eu,e0,alamb,esm,fyd,astotal)
        ff=f;
        prod = fi * ff;
#        print("qi =",qi,"fi =",fi)
#        print("qi =",qf,"fi =",ff)
#
#
#  O intervalo solução foi definido
#  A linha neutra qsi fica entre [qi,qf]
#
#
#
#  O intervalo solução foi definido
#  A linha neutra qsi fica entre [qi,qf]
#  
#  Processo iterativo da bissecante
    fk = 1;
    while np.abs(fk) > 0.0001:
#          
        qk = (qi * ff - qf * fi) / (ff - fi);
        qsi = qk;
        f,nrd,mxrd,myrd=funcao(nsd,qsi,alfa,nc1,xc,yc,alfac,sigmacd,ns,xs,ys,eu,e0,alamb,esm,fyd,astotal)
        fk=f
        prod = fk * fi;
        if prod >= 0:
#
            qi = qk;
            fi = fk;
#                
        else:
#               
            qf = qk;
            ff = fk;
#
#
#  Convergência alcançada
#  qk é a raiz da função f(qsi) dada na equação (2.5.11) do Volume 3 de Curso de Concreto Armado
#
#    print(f,nrd,mxrd,myrd)
    theta[k]=alfa
    nid[k]=nrd
    mixd[k]=mxrd
    miyd[k]=myrd
#
#  Diagrama de interação da seção
#    
plt.figure(figsize=(8.5, 6))
plt.plot(mixd,miyd)
plt.grid()
plt.xlabel("Mxrd")
plt.ylabel("Myrd")

plt.show()    
    
    

  


