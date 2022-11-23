function plotbj(name1, name2)

x1 = load(name1);
x2 = load(name2);

plot(x1(:,1),abs(complex(x1(:,end-1),x1(:,end))),x2(:,1),abs(complex(x2(:,end-1),x2(:,end))))

end

