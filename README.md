# ESP32AlgebraFilters

## Matrix and Vector Library for ESP32 and Embedded Systems

This C++ library for PlatformIO is designed for embedded applications. Originally developed for the ESP32, it is optimized for efficiency without hardware vectorization, making it well-suited for CPU-based calculations. Flexible and easy to use, it meets the needs of robotics projects by providing lightweight support for basic matrix operations and state estimation.

## Key Features
- **Type Agnostic Vectors**: The vector class supports any data type, allowing for calculations on various types when applicable.
- **Optimized Matrix Operations**: Efficient handling of symmetric matrices and decomposition methods (e.g., Cholesky decomposition) for performance improvements.
- **Modular and Extensible**: Built with modularity in mind, making it easy to integrate with other projects or expand with additional functionality.

## Master branch
This branch holds the latest stable release of the library. It is based on the most recent fully tested branch [filter-dev](https://github.com/robinAZERTY/ESP32AlgebraFilters/tree/filter-dev), including a specific version of linearAlgebra, not obviously the latest one.


## Usage
This library is ideal for anyone learning how filters work in C++ on embedded platforms. The current release includes foundational matrix and vector classes, a Kalman filter, and a UKF.

## Roadmap
- validate the first version of the library by passing tests

## Developer's Note
I'm a student working on this library out of curiosity and a desire to learn about matrix computations and state estimators. Through this project, I aim to deepen my understanding of these concepts, particularly in the context of embedded systems.

---
_Note: This README was translated from French to English with the assistance of ChatGPT._
