
Y = DCM.xY.y{1};
try X = Hc{1}; end
try l = DCM.xY.name; end

n = 0;

for i = 1:size(Y,2)
    for j = 1:size(Y,3)
        n = n + 1;
        
        subplot(size(Y,2),size(Y,3),n), plot(squeeze(real(Y(:,i,j))),'b'); hold on;
                              try       plot(squeeze(real(X(:,i,j))),'r');  end
                              try l;
                                  if i == 1 || j == 1 || i == j
                                  t = ([l{i} ' * ' l{j}]);
                                  title(t,'fontsize',16);
                                  end
                              end
                              if i == size(Y,2) && j == size(Y,3)
                                  xlabel('Frequency (Hz)','fontsize',16);
                              end

        
        if i == j;
            plot(squeeze(real(Y(:,i,i))),'g'); hold on;
       try  plot(squeeze(real(X(:,i,i))),'r'); end
        end
    end
end
