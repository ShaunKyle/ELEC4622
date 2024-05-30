# Style guide

Uh... Don't treat any single voice as gospel, I guess?

## C style

I skim read this [C programming practices document](https://github.com/mcinglis/c-style) 
by Malcolm Inglis, and I thought it was written with a good set of guiding 
principles. I'm liking the 80 char line limit.

A post on Chris Wellon's (null program) blog, 
[My personal C coding style as of late 2023](https://nullprogram.com/blog/2023/10/08/), 
covers the reasoning behind some personal style choices he makes. I think the 
practices related to structures (always typedef, prefer struct returns instead 
of output parameters) are worth considering.

## C++ style

idk. A good starting point for modern C++ (C++11 and newer) would be the 
[C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines).

What is modern C++ anyway? As far as I can tell, it means making less C-style 
bugs by writing less C. The article 
[Welcome back to C++](https://learn.microsoft.com/en-us/cpp/cpp/welcome-back-to-cpp-modern-cpp) 
from Microsoft Learn gives a nice overview of modern C++ features. Here's my summary:

Follow the principle of RAII:

- Prevent memory leaks by using smart pointer types such as `std::unique_ptr` 
to handle allocation and deletion of memory (as opposed to manually handling 
memory with `new` and `delete` or `malloc()` and `free()`).
- Reduce unnecessary memory copies by using move semantics.
- Prefer to use standard library containers such as `std::vector` that follow 
the principle of RAII, rather than using C-style arrays or implementing custom 
data structures.

Use less C:

- Eliminate C-style strings (`\0` termination) by using `std::string` instead.
- Eliminate C-style macros by using `constexpr` variables instead.
- Eliminate type-unsafe function pointers by using inline lambda expressions.

There are also some nice quality-of-life language features, like range-based 
for loops.

## When should I choose to use C or C++?

... dunno. [Let's ask Linus Torvalds](http://harmful.cat-v.org/software/c++/linus).

A hardcore C programmer will probably tell you that C++ is a convoluted mess 
of a programming language filled with footguns, and they'd be right.

If you want to use `<insert language feature>`, use C++. Just keep in mind that 
the language continues to grow, and other developers might not use the same 
subset of features that you use. If there is another programming language which 
matches the set of features you'd like, consider using that.

If you want to write in a systems programming language that won't change much 
over time, use C. Just keep in mind that there are common classes of bugs that 
are hard to avoid, and language features that you might take for granted won't 
be present. There are other options for systems programming languages, but C is 
still the lingua franca.

Don't attempt to use C++ as "C, but with `<insert language feature>`". That 
just gives you the worst of both worlds.

## When should I give up and use a memory-managed language instead?

Depends on how much of a masochist you are.

## Is there a viable alternative to C and C++?

Forth.
