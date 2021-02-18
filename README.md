# amath
* c++ math library ( currently on work )
# Documentation

## Points
* You can create n-dimensional points.
>Points p1({1, 2, 3, 4});
>
>cout << "Point is :  " << p1;
>
>>**Point is :  ( 1, 2, 3, 4, )**

* You can change the value of an axis after initialization.
>Points p2({1, 2, 3, 4});
>
>p2.setAxis(1,4);
>
>cout << "Point is :  " << p2;
>
>>**Point is :  ( 1, 4, 3, 4, )**

* You can get the value of an axis.

>cout << "Field is :  " << p2.getAxis(1);
>
>>**Field is :   4**
>>

* You can sum, subtract points even if they are **not the same size.**
>cout << "Point is :  " << p1 + p2 << endl;
>
>cout << "Point is :  " << p1 - p2;
>
>>**Point is :  ( 2, 6, 6, 8, )**
>>
>>**Point is :  ( 0, -2, 0, 0, )**
>>

* You can increment and decrement points.
>cout << "Point is :  " << p1++;
>


>>**Point is :  ( 2, 3, 4, 5, )**

* You can get the distance between points. 

>Point p3({1, 2}); Point p4({1});
>
>cout << "Distance is  :  " << p3.distance(p4);
>
>>**Distance is  :  2**
>>
>
* You don't have to create points before hand to make them work.
>Point p3({1, 2}); 
>
>cout << "Distance is  :  " << p3.distance(Point({1}));
>
>>**Distance is  :  2**

* You can output and input points directly as shown before.
## Matrix
* You can create matrices with desired size and shape, you can modify shape of matrices, take transpose of matrices.
* You can subtract, sum, multiply matrices with each other.
* You can output matrices directly to the screen in compliance with specified shape.
* .size() to get element count.
>Matrix m1(2,3,{
1, 2,
3, 4,
5, 6
});
>Matrix m2 = m1.transpose();
>
>Matrix m3 = m1 * m2;
>
>Matrix m4 = m1 + m1;
>
>Matrix m5 = m4 * 4;
>
>Matrix m6 = m5.reShape(6, 1);
>
>cout << m1.size()
>
>cout << m1 << endl;
>
## Poly
* You can derivate and integrate polynomials **definitely** and **numerically**.
* You can evaluate them at a certain coordinate and this will return an *Point*.
* You can change a coefficient or order of a term after initialization. (*currently this will cause problems with summing and subtracting with other polynomials*)
* You can get the term count.
* You can output polynomials directly with polynomial display.

## Shape(currently on work)
## Functions(currently on work)
