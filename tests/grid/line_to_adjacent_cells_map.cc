// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <iterator>
#include <set>

#include "../tests.h"


using namespace dealii;


void
check(Triangulation<3> &tria)
{
  const auto vertex_to_cells = GridTools::vertex_to_cell_map(tria);

  for (const auto &cell : tria.active_cell_iterators())
    for (const unsigned int line : cell->line_indices())
      Assert(cell->get_cells_adjacent_to_line(line).empty(),
             ExcInternalError());

  for (unsigned int repetition = 0; repetition < 2; ++repetition)
    {
      tria.compute_line_to_adjacent_cells_map();

      for (const auto &cell : tria.active_cell_iterators())
        for (const unsigned int line : cell->line_indices())
          {
            const auto  edge = cell->line(line);
            const auto &a    = vertex_to_cells[edge->vertex_index(0)];
            const auto &b    = vertex_to_cells[edge->vertex_index(1)];
            std::set<Triangulation<3>::active_cell_iterator> expected;
            std::set_intersection(a.begin(),
                                  a.end(),
                                  b.begin(),
                                  b.end(),
                                  std::inserter(expected, expected.end()));

            Assert(cell->get_cells_adjacent_to_line(line) == expected,
                   ExcInternalError());
          }
    }
}


int
main()
{
  initlog();

  Triangulation<3> hex;
  GridGenerator::subdivided_hyper_cube(hex, 2);
  check(hex);
  hex.begin_active()->set_refine_flag();
  hex.execute_coarsening_and_refinement();
  check(hex);

  Triangulation<3> tet;
  GridGenerator::subdivided_hyper_cube_with_simplices(tet, 2);
  check(tet);

  deallog << "OK" << std::endl;
}
