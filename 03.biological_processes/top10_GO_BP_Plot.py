import pandas as pd
import matplotlib.pyplot as plt

def plot_top10(ax, file, fill_color, title_text):
    df = pd.read_csv(file)
    df = df.sort_values("Binom_Fold_Enrichment", ascending=False).head(10)
    df = df.sort_values("Binom_Fold_Enrichment", ascending=True)
    ax.barh(df["name"], df["Binom_Fold_Enrichment"], color=fill_color)
    ax.set_xlabel("Fold Enrichment")
    ax.set_title(title_text)

fig, axes = plt.subplots(3, 1, figsize=(12, 18))
plot_top10(axes[0], "results/all_human_BP_filtered.csv","#6a8caf", "Human All OCRs (hg38)")
plot_top10(axes[1], "results/all_mouse_BP_filtered.csv","#a8b87a", "Mouse All OCRs (mm10)")
plot_top10(axes[2], "results/human_shared_BP_filtered.csv","#a47db8", "Shared OCRs (hg38 coordinates)")
plt.tight_layout(pad=3.0)
plt.savefig("results/plots_combined.png", dpi=150, bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(figsize=(12, 6))
plot_top10(ax, "results/human_specific_BP_filtered.csv", "#c87941", "Human Specific OCRs (hg38)")
plt.tight_layout()
plt.savefig("results/plots_human_specific.png", dpi=150, bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(figsize=(12, 6))
plot_top10(ax, "results/mouse_specific_BP_filtered.csv", "#c47a8a", "Mouse Specific OCRs (mm10)")
plt.tight_layout()
plt.savefig("results/plots_mouse_specific.png", dpi=150, bbox_inches="tight")
plt.close()