nodes = [1 0;
         1 1;
         0 1;
         0 0;
         .5 .5;
         1 .5;
         .5 1;
         0 .5;
         .5 0;
         .75 .25;
         .75 .75;
         .25 .75;
         .25 .25] .*10

elem = [1 4 5 9 13 10;
        1 2 5 6 11 10;
        2 3 5 7 12 11;
        3 4 5 8 13 12]
edge = [4 9 1;
        1 6 2;
        2 7 3;
        3 8 4]
bc = [2 -1 4 -1]'
f(x) = [[1e6 0]',[0 0]']
