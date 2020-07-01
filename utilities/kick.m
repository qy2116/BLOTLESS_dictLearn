function pos = kick(SP)
n = length(SP);
pos = [];
for ind = 1:n
    if isempty(SP{ind})
        pos = [pos ind];
    end
end
end