clc
clear
nomvars=['Numero de Mezclas a cuantificar';'Número de sustancias a cuantificar';'Concentración máxima en espectros de calibración ';'Número de puntos de calibración';'Volumen de alícuota (mL)';'Volumen de aforo (mL)';'Ingrese longitudes de onda del espectrofotometro';'Ingrese número de longitudes de onda a seleccionar';'Ingrese repeticiones en los puntos de calibración'];
txtvars=string(zeros(9,1));
txtvars=x_mdialog('Ingrese los siguientes datos',nomvars,txtvars);
valvars=evstr(txtvars);


//Mezclas a Cuantificar
M=valvars(1,1)+1;
//Sustancias a cuantificar
Ns=valvars(2,1);
//Concentración de la mayor concentración
Cmax=valvars(3,1);
//Número de mezclas de calibración
Mcal=valvars(4,1)+1;
//Volumen de la alícuota en mL
Al=valvars(5,1);
//Volumen de aforo en mL
Af=valvars(6,1);
//Número de longitudes de onda medidas por el espectrofotómetro
F=valvars(7,1);
//Número de longitudes de onda seleccionadas
NL=valvars(8,1);
//Repeticiones en los puntos de la curva de calibración
T=valvars(9,1);
//Matriz de datos de calibraciòn: "Cal"
Cal=[zeros(F,Mcal)];
titulos=['Ingrese absorbancia de las mezclas de calibración';'Longitudes de onda en la primera columna'];
Cal=x_matrix(titulos,Cal);
//Matriz de datos de absorbancia de mezclas a cuantificar
X=[zeros(F,M)];
titulos=['Ingrese absorbancia de las mezclas a cuantificar';'Longitudes de onda en la primera columna'];
X=x_matrix(titulos,X);
//Vector de longitudes de onda seleccionadas para calibración y cuantificación
Y=[zeros(NL,1)];
titulos=['Ingrese Longitudes de onda a seleccionar'];
Y=x_matrix(titulos,Y);
//Matriz de concentraciones conocidas de cada mezcla de calibración
Cc=[zeros(Ns,Mcal-1)];
titulos=['Ingrese las concentraciones de cada sustancia en cada mezcla de calibración (ug/g)';'Cada fila es una sustancia y cada columna una mezcla diferente']
Cc=x_matrix(titulos,Cc)/Cmax;
//Matriz de las mezclas a cuantificar con valores de absorbancia seleccionados
Z=[zeros(NL,M)]
for i=1:NL
    for j=1:F
    if X(j,1)==Y(i) 
        Z(i,:)=X(j,:) 
    end
    end
   end
//Matriz de las mezclas de calibración con valores de absorbancia seleccionados
R=[zeros(NL,Mcal)]
for i=1:NL
    for j=1:F
    if Cal(j,1)==Y(i) 
        R(i,:)=Cal(j,:)
    end
    end
   end
//Matriz de espectros de las sustancias en mezcla "S"
S=(inv(Cc*Cc')*Cc*R(1:NL,2:Mcal)')';
Sf=inv(S'*S)*S';
Cd=[zeros(M-1,Ns)];
for i=1:M-1
for j=1:Ns
    if j==1
Cd(i,j)=Sf(j,:)*Z(:,i+1)*Cmax*Af/Al;
        if Cd(i,j)<0
            Cd(i,j)=0
            end
else
    if j==2
    Cd(i,j)=Sf(j,:)*Z(:,i+1)*Cmax*Af/Al;
            if Cd(i,j)<0
            Cd(i,j)=0
        end
    else 
            Cd(i,j)=Sf(j,:)*Z(:,i+1)*Cmax*Af/Al;
            if Cd(i,j)<0
            Cd(i,j)=0
        end
        end
end
end
end
i=1
j=1
editvar Cd
//CALCULO DE FIGURAS DE MÉRITO
for i=1:NL
k=2
for j=1:(Mcal-1)/T
        Desv(i,j)=stdev(R(i,k:k+2))
        Rp(i,j)=mean(R(i,k:k+2))
        k=k+T
end
end

for j=1:(Mcal-1)/T
    De(j)=norm(Desv(:,j))
    NRp(j)=norm(Rp(:,j))
    SN(j)= NRp(j)/De(j)
    end

i=1
j=1
for i=1:Ns
    //Cálculo de sensibilidad a ada sustancia (absorbancia/masa)
    SEN(i)=1/(norm(Sf(i,:))*Cmax);
    //Cálculo de la selectividad 
    SEL(i)=SEN(i)*Cmax/(norm(S(:,i)));
    //Cálculo del límite de detecciòn
    LOD(i)=3*Af*mean(De)/(SEN(i)*Al);
end

disp('Selectividad, Sensibilidad, Relación Señal/Ruido, Límite de detección')
LOD 
SEN 
SEL
//I=Mcal
//J=NL
//K=Ns
//R=Z
//end
Rm=inv(R(1:NL,2:Mcal)'*R(1:NL,2:Mcal))*R(1:NL,2:Mcal)'; //(24x540); R (540x25)
b=Rm'*Cc';//(540x3)
I=eye(NL,Mcal-1)
varc=R'*(I-R(1:NL,2:Mcal)+R(1:NL,2:Mcal))*R/(NL-Ns)
//varr=
for i=1:M
    h(i)=R(:,i)'*Rm'
end
