# IKproto

We provided two algorithms in this repository dedicated to the robotic arm given below.

![335736331_242841668070773_6298975862098763900_n](https://user-images.githubusercontent.com/83645103/226197082-3bd30e9d-580a-433e-be26-4f288e85383c.png)

First "dh" which one computes forward kinematics transfomartion by Denavid-Harthenberg parameters symbollically and numerically.

Second, "main" indicates an inverse kinematics algorithm to compute rotation angles of all axes in an original way by rotation and translation matrices and analytical calculations.
Inputs of second algorithm are:

  xb, yb, zb, rXb, rYb, rZb = position and rotation of the robot's 0 coordinate system in reference to the basic simulation system
  
  x4, y4, z4, rX4, rY4, rZ4 = position and rotation of the robot's tool coordinate system in reference to the basic simulation system
  
For a more digital explanation in Polish, see the pdf below:

[IK_visible_model_v3.0.pdf](https://github.com/RadoslawDebinski/IKproto/files/11011911/IK_visible_model_v3.0.pdf)
