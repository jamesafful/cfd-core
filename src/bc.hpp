
#pragma once
enum class BCType { Inflow, Outflow, Wall, Periodic };
struct BC { BCType left{BCType::Periodic}, right{BCType::Periodic}; };
