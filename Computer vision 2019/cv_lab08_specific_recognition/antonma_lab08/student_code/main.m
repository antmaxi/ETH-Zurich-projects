load("dataset.mat");
i = 1;
j = 2;
display_flag = 1;
X1 = objects(i).X;
X2 = objects(j).X;

matchingCost = shape_matching(X1,X2,display_flag);
