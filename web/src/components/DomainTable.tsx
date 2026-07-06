import { useMemo, useState } from "react";
import {
  createColumnHelper,
  flexRender,
  getCoreRowModel,
  getFilteredRowModel,
  getSortedRowModel,
  useReactTable,
  type SortingState,
} from "@tanstack/react-table";
import type { DomainCall } from "../lib/types";

const helper = createColumnHelper<DomainCall>();

const COLUMNS = [
  helper.accessor("status", { header: "Status", sortingFn: "alphanumeric" }),
  helper.accessor("gene", { header: "Gene", sortingFn: "alphanumeric" }),
  helper.accessor("name", { header: "Domain", sortingFn: "alphanumeric" }),
  helper.accessor("accession", {
    header: "Accession",
    cell: (info) => <code>{info.getValue()}</code>,
    sortingFn: "alphanumeric",
  }),
  helper.accessor("type", { header: "Type", sortingFn: "alphanumeric" }),
  helper.accessor("start", {
    header: "Range",
    cell: (info) => `${info.getValue()}–${info.row.original.end}`,
    sortingFn: "basic",
    // filter on the rendered "start–end" string
    filterFn: (row, _colId, filterValue: string) => {
      const val = `${row.original.start}–${row.original.end}`;
      return val.includes(filterValue);
    },
  }),
];

// Default sort: gene asc, then start asc
const DEFAULT_SORTING: SortingState = [
  { id: "gene", desc: false },
  { id: "start", desc: false },
];

interface Props {
  domains: DomainCall[];
}

export function DomainTable({ domains }: Props) {
  const [sorting, setSorting] = useState<SortingState>(DEFAULT_SORTING);
  const [globalFilter, setGlobalFilter] = useState("");

  const data = useMemo(() => domains, [domains]);

  const table = useReactTable({
    data,
    columns: COLUMNS,
    state: { sorting, globalFilter },
    onSortingChange: setSorting,
    onGlobalFilterChange: setGlobalFilter,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
    globalFilterFn: "includesString",
  });

  const sortIcon = (colId: string) => {
    const col = table.getColumn(colId);
    if (!col?.getIsSorted()) return " ↕";
    return col.getIsSorted() === "asc" ? " ↑" : " ↓";
  };

  return (
    <div className="domain-table-wrapper">
      <div className="domain-table-toolbar">
        <input
          className="domain-table-filter"
          type="search"
          placeholder="Filter domains…"
          value={globalFilter}
          onChange={(e) => setGlobalFilter(e.target.value)}
          aria-label="Filter domain table"
        />
        <span className="domain-table-count">
          {table.getFilteredRowModel().rows.length} /{" "}
          {domains.length} rows
        </span>
      </div>
      <table className="domain-table">
        <caption className="domain-table-caption">
          Raw domain calls from InterPro/Pfam — the diagram merges overlapping
          entries of the same region and only renders{" "}
          <em>domain</em>, <em>repeat</em>, and <em>conserved_site</em> types.
          Click a column header to sort.
        </caption>
        <thead>
          {table.getHeaderGroups().map((hg) => (
            <tr key={hg.id}>
              {hg.headers.map((header) => (
                <th
                  key={header.id}
                  onClick={header.column.getToggleSortingHandler()}
                  style={{ cursor: header.column.getCanSort() ? "pointer" : "default" }}
                  aria-sort={
                    header.column.getIsSorted() === "asc"
                      ? "ascending"
                      : header.column.getIsSorted() === "desc"
                        ? "descending"
                        : "none"
                  }
                >
                  {flexRender(header.column.columnDef.header, header.getContext())}
                  {header.column.getCanSort() && (
                    <span className="sort-icon" aria-hidden="true">
                      {sortIcon(header.column.id)}
                    </span>
                  )}
                </th>
              ))}
            </tr>
          ))}
        </thead>
        <tbody>
          {table.getRowModel().rows.map((row) => (
            <tr
              key={row.id}
              className={`status-${row.original.status.toLowerCase()}`}
            >
              {row.getVisibleCells().map((cell) => (
                <td key={cell.id}>
                  {flexRender(cell.column.columnDef.cell, cell.getContext())}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
