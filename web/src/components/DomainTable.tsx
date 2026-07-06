import { useMemo, useState } from "react";
import {
  createColumnHelper,
  flexRender,
  getCoreRowModel,
  getFacetedRowModel,
  getFacetedUniqueValues,
  getFilteredRowModel,
  getSortedRowModel,
  useReactTable,
  type Column,
  type ColumnFiltersState,
  type SortingState,
} from "@tanstack/react-table";
import type { DomainCall } from "../lib/types";

// Columns that should render a <select> filter instead of a text input.
const SELECT_COLUMNS = new Set(["status", "gene", "type"]);

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
    filterFn: (row, _colId, filterValue: string) =>
      `${row.original.start}–${row.original.end}`.includes(filterValue),
  }),
];

// Default sort: gene asc → start asc
const DEFAULT_SORTING: SortingState = [
  { id: "gene", desc: false },
  { id: "start", desc: false },
];

// Numeric superscript badges for sort priority (up to 3 active sorts)
const PRIORITY_BADGE = ["①", "②", "③"];

function SortHeader({ column }: { column: Column<DomainCall, unknown> }) {
  const sorted = column.getIsSorted();
  const sortIndex = column.getSortIndex();
  return (
    <span className="sort-header-content">
      {sorted && sortIndex >= 0 && sortIndex < PRIORITY_BADGE.length && (
        <span className="sort-priority" aria-hidden="true">
          {PRIORITY_BADGE[sortIndex]}
        </span>
      )}
      <span className="sort-arrow" aria-hidden="true">
        {sorted === "asc" ? "↑" : sorted === "desc" ? "↓" : "↕"}
      </span>
    </span>
  );
}

/** Dropdown populated from distinct values currently in the (pre-filter) data. */
function SelectFilter({ column }: { column: Column<DomainCall, unknown> }) {
  const unique = Array.from(column.getFacetedUniqueValues().keys()).sort();
  const value = (column.getFilterValue() ?? "") as string;
  return (
    <select
      className="col-filter col-filter-select"
      value={value}
      onChange={(e) => column.setFilterValue(e.target.value || undefined)}
      aria-label={`Filter ${column.id}`}
    >
      <option value="">All</option>
      {unique.map((v) => (
        <option key={v} value={v}>
          {v}
        </option>
      ))}
    </select>
  );
}

/** Text input for free-text columns. */
function TextFilter({ column }: { column: Column<DomainCall, unknown> }) {
  const value = (column.getFilterValue() ?? "") as string;
  return (
    <input
      className="col-filter col-filter-text"
      type="search"
      placeholder="filter…"
      value={value}
      onChange={(e) => column.setFilterValue(e.target.value || undefined)}
      aria-label={`Filter ${column.id}`}
    />
  );
}

interface Props {
  domains: DomainCall[];
}

export function DomainTable({ domains }: Props) {
  const [sorting, setSorting] = useState<SortingState>(DEFAULT_SORTING);
  const [columnFilters, setColumnFilters] = useState<ColumnFiltersState>([]);

  const data = useMemo(() => domains, [domains]);

  const table = useReactTable({
    data,
    columns: COLUMNS,
    state: { sorting, columnFilters },
    onSortingChange: setSorting,
    onColumnFiltersChange: setColumnFilters,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
    getFacetedRowModel: getFacetedRowModel(),
    getFacetedUniqueValues: getFacetedUniqueValues(),
    enableMultiSort: true,
    maxMultiSortColCount: 3,
  });

  const filtered = table.getFilteredRowModel().rows.length;
  const hasFilters = columnFilters.length > 0;

  return (
    <div className="domain-table-wrapper">
      <div className="domain-table-toolbar">
        <span className="domain-table-count">
          {filtered !== domains.length
            ? `${filtered} of ${domains.length} rows`
            : `${domains.length} rows`}
        </span>
        {hasFilters && (
          <button
            type="button"
            className="clear-filters-btn"
            onClick={() => setColumnFilters([])}
          >
            Clear filters
          </button>
        )}
        <span className="domain-table-hint">
          Click header to sort · Shift+click for multi-column sort
        </span>
      </div>
      <table className="domain-table">
        <caption className="domain-table-caption">
          Raw domain calls from InterPro/Pfam — the diagram merges overlapping
          entries of the same region and only renders{" "}
          <em>domain</em>, <em>repeat</em>, and <em>conserved_site</em> types.
        </caption>
        <thead>
          {/* Sort row */}
          {table.getHeaderGroups().map((hg) => (
            <tr key={hg.id}>
              {hg.headers.map((header) => (
                <th
                  key={header.id}
                  onClick={header.column.getToggleSortingHandler()}
                  style={{ cursor: "pointer" }}
                  aria-sort={
                    header.column.getIsSorted() === "asc"
                      ? "ascending"
                      : header.column.getIsSorted() === "desc"
                        ? "descending"
                        : "none"
                  }
                >
                  <span className="th-label">
                    {flexRender(header.column.columnDef.header, header.getContext())}
                    <SortHeader column={header.column} />
                  </span>
                </th>
              ))}
            </tr>
          ))}
          {/* Per-column filter row */}
          <tr className="filter-row">
            {table.getLeafHeaders().map((header) => (
              <td key={header.id} className="filter-cell">
                {SELECT_COLUMNS.has(header.column.id) ? (
                  <SelectFilter column={header.column} />
                ) : (
                  <TextFilter column={header.column} />
                )}
              </td>
            ))}
          </tr>
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
