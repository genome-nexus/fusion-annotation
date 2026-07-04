import { EXAMPLES } from "../lib/examples";
import type { AnnotateParams } from "../lib/types";
import { DEFAULT_PARAMS } from "../lib/defaultParams";

interface Props {
  onSelect: (params: AnnotateParams) => void;
  disabled: boolean;
}

/** Row of buttons that fill in and immediately run a known-good fusion
 * lookup — a quick way for a first-time visitor to see real output without
 * having to already know a valid gene/exon combination. */
export function ExampleFusions({ onSelect, disabled }: Props) {
  return (
    <div className="examples-row">
      <span className="examples-label">Try an example:</span>
      {EXAMPLES.map((ex) => (
        <button
          key={ex.label}
          type="button"
          className="example-button"
          title={ex.description}
          disabled={disabled}
          onClick={() => onSelect({ ...DEFAULT_PARAMS, ...ex.params })}
        >
          {ex.label}
        </button>
      ))}
    </div>
  );
}
