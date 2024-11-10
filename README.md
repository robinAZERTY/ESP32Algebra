# ESP32AlgebraFilters

## Matrix and Vector Library for ESP32 and Embedded Systems

This repository contains a custom C++ matrix and vector library tailored for use in embedded systems, particularly for robotic applications on platforms like the ESP32. Originally developed to support Kalman filtering in drone and mobile robotics projects, this library is designed for flexibility, efficiency, and ease of use. It includes optimized implementations to minimize computational redundancy, especially useful for resource-constrained environments.

## Key Features
- **Type Agnostic Vectors**: The vector class supports any data type, allowing for calculations on various types when applicable.
- **Optimized Matrix Operations**: Efficient handling of symmetric matrices and decomposition methods (e.g., Cholesky decomposition) for performance improvements.
- **Modular and Extensible**: Built with modularity in mind, making it easy to integrate with other projects or expand with additional functionality.

## Usage
This library is ideal for anyone implementing filters or control algorithms in C++ on embedded platforms. The current release includes foundational matrix and vector classes, and a Kalman filter class will be added soon.

## Roadmap
- Add Kalman filter class for both linear (KF) and nonlinear (UKF) applications.
- Implement unit tests for matrix and vector operations.
- Continue optimizing code for embedded use cases.

## Developper's Note
I'm a student working on this library out of curiosity and a desire to learn about matrix computations and state estimators. Through this project, I aim to deepen my understanding of these concepts, particularly in the context of embedded systems.
