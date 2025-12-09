
#pragma once
enum class BCType { Transmissive, Periodic };
struct BC1D { BCType left{BCType::Transmissive}, right{BCType::Transmissive}; };
