function count = crossing_rate(v,soglia)

count=0;
for k=2:length(v)
    if v(k) >= soglia && v(k-1,1) < soglia
        count=count+1;
    end
end