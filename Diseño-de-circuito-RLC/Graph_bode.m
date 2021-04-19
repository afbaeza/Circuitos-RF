clear all
close all
clc

% % % % % % % % Lee los datos obtenidos de la simulación en formato .RAW 
% % % % % % % % generado por Ltspice y los convierte en arreglos para
% % % % % % % % realizar el tratamiento que se le desee dar.



% Lee el archivo .RAW generado por Ltspice y los guado en Datos
FILENAME = 'Design_RLC.raw';


Datos = LTspice2Matlab(FILENAME);

% Almacena en f el vector que contiene los datos de la frecuencia.
f = Datos.freq_vect(1,:);


% Almacena en Data_mag_phase los valores de magnitud y fase del diagrama de
% bode obtenido por la simulación 
for i = 1:Datos.num_steps
    Data_mag_phase(i,:) = Datos.variable_mat(3,:,i);
end

% Se comvierten todos los valores de magnitud a decibeles y los de la fase
% a grados
[m, n] = size(Data_mag_phase);
for i = 1:m
    for j = 1:n
        Bode_mag(i,j) = abs(Data_mag_phase(i,j));
        Bode_mag(i,j) = 20*log10(Bode_mag(i,j));
        Bode_phase(i,j) = angle(Data_mag_phase(i,j));
        Bode_phase(i,j) = Bode_phase(i,j)*(180/pi);
    end
end

%%
% tener precaución esta parte del algoritmo puede tardarse un mas
% tiempo de lo normal dependiendo de la cantidad de simulaciones
% realizadas 
close all
% % % % % % % % % % % % % % %    Graficas   % % % % % % % % % % % %
% % % % % % % % % % % % % % %               % % % % % % % % % % % % 
% % % % % % % % % % % % % % %               % % % % % % % % % % % % 
% 1k
f0 = 1.49968e8;
BW = 3.60433e6;
fl =   f0 - (BW/2);
fh =   f0 + (BW/2);

subplot(2,1,1)

semilogx(f, Bode_mag(i,:)+20,"LineWidth", 2);
hold on
grid on

title('Diagrama de Bode de simulación')
ylabel('Magnitud [db]')
xlabel('Frecuencia [Hz]')
xline(f0,'--k')
% xline(fl,'--k')
% xline(fh,'--k')
ylim([-100 0])
xlim([0 1e9])

text(f0,-10,'\leftarrow f_0 = 149.96899 [MHz] ')

% legend( "fc < 150 [MHz]")
% legend("boxoff")

hold off 

subplot(2,1,2)
semilogx(f, Bode_phase(i,:),"LineWidth", 2);
hold on
grid on
xlabel('Frecuencia [Hz]')
ylabel('Fase [Grados]')
hold off 
%%

% % % Guarda las Graficas en formato .fig y .eps

savefig("Bodesim.fig")
saveas(gcf, "Bodesim.eps", "epsc")


disp('Finalizado')