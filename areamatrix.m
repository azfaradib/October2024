[num,txt,raw] = xlsread('Data mean std.xlsx');

for i=1:(size(num,1))
   average(i)=num(i,2);
    stddev(i)=num(i,3);
end


for i=1:1
    a=average(i);
    b=stddev(i);
    for j=1:1
      c= average(j);
      d=stddev(j);
output=area_between(a,b,c,d);
   result(i,j)=output;
    end
end



