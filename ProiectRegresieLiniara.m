x1_id = id.X{1};
x2_id = id.X{2};
y_id = id.Y;

x1_val = val.X{1};
x2_val = val.X{2};
y_val = val.Y;

Lx1_id = length(x1_id);
Lx2_id = length(x2_id);

Lx1_val = length(x1_val);
Lx2_val = length(x2_val);

figure(1)
mesh(x1_id, x2_id, y_id)
title('Datele de identificare')

figure(2)
mesh(x1_val, x2_val, y_val)
title('Datele de validare')

m = 30; 
MSE_id = zeros(1, m);
MSE_val = zeros(1, m);

for grad_polinom = 1:m
    phi_id = phi_f(grad_polinom, x1_id, x2_id);
   
    y_idfinal = reshape(y_id, [], 1);
    teta = phi_id \ y_idfinal;
  
    y_id_aprox = phi_id * teta;
    Y_Aproximat_id = reshape(y_id_aprox, Lx1_id, Lx2_id);

    valoareMSE_id = (y_id_aprox - y_idfinal).^2;
    mse_id = mean(valoareMSE_id);  
    MSE_id(grad_polinom) = mse_id;
    
    phi_val = phi_f(grad_polinom, x1_val, x2_val);

    y_valfinal = reshape(y_val, [], 1);
    y_val_aprox = phi_val * teta;
    Y_Aproximat_val = reshape(y_val_aprox, Lx1_val, Lx2_val);

    valoareMSE_val = (y_val_aprox - y_valfinal).^2;
    mse_val = mean(valoareMSE_val);  
    MSE_val(grad_polinom) = mse_val;
    
end

MSE_val_min = min(MSE_val);
for i = 1:m
    if MSE_val_min == MSE_val(i)
        z = i;
    end
end

MSE_id_min = min(MSE_id);
for i = 1:m
    if MSE_id_min == MSE_id(i)
        k = i;
    end
end

%Aproximare identificare
for i=1:k
    phi = phi_f(k, x1_id, x2_id);
    y_idfinal = reshape(y_id, [], 1);
    teta = phi \ y_idfinal;
    y_id_aprox = phi * teta;
    Y_Aproximat_id = reshape(y_id_aprox, Lx1_id, Lx2_id);

end
figure(3)
mesh(x1_id, x2_id, Y_Aproximat_id)
title('Aproximare date de identificare')
xlabel('X1')
ylabel('X2')
zlabel('Y')

%Aproximare validare
for i=1:z
    phi = phi_f(z, x1_id, x2_id);
    y_idfinal = reshape(y_id, [], 1);
    teta = phi \ y_idfinal;

    phi_val = phi_f(z, x1_val, x2_val);

    y_val_aprox = phi_val * teta;
    Y_Aproximat_val = reshape(y_val_aprox, Lx1_val, Lx2_val);
end

figure(4)
mesh(x1_val, x2_val, Y_Aproximat_val)
title('Aproximare date de validare');
xlabel('X1')
ylabel('X2')
zlabel('Y')

fprintf('Valoare MSE_id_min: %f, la gradul: %d', MSE_id_min, k)
fprintf('\nValoare MSE_val_min: %f, la gradul: %d\n', MSE_val_min, z)

figure(5)
plot(MSE_id, '-o'), grid;
title('MSE identificare')


figure(6)
plot(MSE_val, '-o'); grid;
title('MSE ')



function phi = phi_f(grad_polinom, x1, x2)
    nr_t = (grad_polinom + 2) * (grad_polinom + 1) / 2; 
    linia = 1;
    phi = zeros(length(x1) * length(x2), nr_t);
    
    for i = 1:length(x1)
        for j = 1:length(x2)
            coloana = 1;
            for a = 0:grad_polinom
                for b = 0:(grad_polinom - a)
                    phi(linia, coloana) = (x1(i)^a) * (x2(j)^b);
                    coloana = coloana + 1;
                end
            end
            linia = linia + 1;
        end
    end
end

