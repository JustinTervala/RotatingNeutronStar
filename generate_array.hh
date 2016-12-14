
/************************************************************/
/* Static Compile-time Table Generation of Grid using Variadic Templates */
/************************************************************/

//Generate sequence 
template<std::size_t... Is> struct seq{};

template<std::size_t N, std::size_t... Is> 
struct gen_seq : gen_seq<N-1, N-1, Is...> {};

template<std::size_t... Is>
struct gen_seq<0, Is...> : seq<Is...>{};

//Generator Function
template<class Generator, std::size_t...Is>
constexpr auto generate_array_helper(Generator g, seq<Is...>)
-> std::array<decltype(g(std::size_t{}, sizeof...(Is))), sizeof...(Is)>
{
    return {{g(Is, sizeof...(Is))...}};
}

template<std::size_t tcount, class Generator>
constexpr auto generate_array(Generator g)
-> decltype(generate_array_helper(g, gen_seq<tcount>{}))
{
    return generate_array_helper(g, gen_seq<tcount>{});
}

/* 
 * USEAGE EXAMPLE
 * 
 * struct point {
 *     double x, y;
 * };
 *
 * #include <iostream>
 * std::ostream& operator<<(std::ostream& o, const point& p) 
 * { return o << p.x << ", " << p.y; }
 *
 * //User-defined generator
 * constexpr point my_generator(std::size_t curr, std::size_t total) {
 *     return {curr*40.0/(total-1), curr*20.0/(total-1)};
 * }
 *
 * int main() {
 *     constexpr auto first_array = generate_array<5>(my_generator);
 *     constexpr auto second_array = generate_array<10>(my_generator);
 *     for (auto p : first_array) {
 *         std::cout << p << '\n';
 *     }
 * }

