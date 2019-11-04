function [myGrid, delX] = create_grid(L,N)
    delX = L/N;
    myGrid = [0, (0 + delX/2):delX:(L - delX/2), L];
end