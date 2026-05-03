"""
Скрипт для чтения и обработки фМРТ данных из датасета Attention
Обновленная версия: добавлена временная фильтрация и генерация иллюстраций для диплома
"""


import os
import glob
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import seaborn as sns
from nilearn import datasets, plotting
from nilearn.input_data import NiftiLabelsMasker



DATA_PATH = r'...\attention\functional'
TR = 2.0  # Время повторения (Repetition Time). 

if not os.path.exists(DATA_PATH):
    print(f"ОШИБКА: Папка {DATA_PATH} не найдена!")
    raise FileNotFoundError(f"Папка не найдена: {DATA_PATH}")
else:
    print(f"✓ Папка найдена: {DATA_PATH}")


print("\n" + "="*70)
print("ШАГ 1: Поиск и загрузка файлов фМРТ")
print("="*70)

patterns_to_try = ['snrfM*.img', 'srfM*.img', 'rfM*.img', 'fM*.img', '*.img']

fmri_files = []
for pattern in patterns_to_try:
    file_pattern = os.path.join(DATA_PATH, pattern)
    fmri_files = sorted(glob.glob(file_pattern))
    if len(fmri_files) > 0:
        print(f"\n✓ Найдено файлов по паттерну '{pattern}': {len(fmri_files)}")
        break

if len(fmri_files) == 0:
    raise FileNotFoundError("Не найдены файлы фМРТ для обработки")


print(f"\nЗагрузка {len(fmri_files)} файлов...")
all_scans = []

for i, file_path in enumerate(fmri_files):
    img = nib.load(file_path)
    data = img.get_fdata()
    if data.ndim == 4 and data.shape[-1] == 1:
        data = data[..., 0]
    all_scans.append(data)

fmri_4d = np.stack(all_scans, axis=-1)
print(f"✓ Все файлы загружены! Итоговая размерность: {fmri_4d.shape}")


img_4d = nib.Nifti1Image(fmri_4d, affine=img.affine)


print("\n" + "="*70)
print("ШАГ 3: Загрузка атласа AAL")
print("="*70)

aal_atlas = datasets.fetch_atlas_aal()
atlas_filename = aal_atlas['maps']



print("\n" + "="*70)
print("ШАГ 4: Извлечение временных рядов (с фильтрацией и стандартизацией)")
print("="*70)


masker = NiftiLabelsMasker(
    labels_img=atlas_filename,
    standardize='zscore', # Z-стандартизация
    detrend=True,         # Удаление линейного тренда
    low_pass=0.1,         # Верхняя граница фильтрации
    high_pass=0.01,        # Нижняя граница фильтрации
    t_r=TR,               # Время повторения
    verbose=1
)

# Для визуализации "До/После" создадим маскер БЕЗ фильтрации
raw_masker = NiftiLabelsMasker(labels_img=atlas_filename)

print("Применяем фильтрацию...")
time_series = masker.fit_transform(img_4d)
raw_time_series = raw_masker.fit_transform(img_4d)

print(f"\n✓ Временные ряды извлечены! Размерность: {time_series.shape}")


variances = np.var(time_series, axis=0)
threshold = 0.01
good_regions = np.where(variances >= threshold)[0]
X = time_series[:, good_regions]
labels = [aal_atlas['labels'][i] for i in good_regions]

output_dir = 'output'
os.makedirs(output_dir, exist_ok=True)
np.save(os.path.join(output_dir, 'fmri_timeseries.npy'), X)

print("\n" + "="*70)
print("ШАГ 8: Генерация иллюстраций для диплома")
print("="*70)

# Рис 4.1: Визуализация атласа
plotting.plot_roi(atlas_filename, title="AAL Atlas Parcellation", 
                  display_mode='ortho', cut_coords=[0, 0, 0])
plt.savefig(os.path.join(output_dir, 'figure_4_1_atlas.png'), dpi=300)
print("✓ Сохранен Рис 4.1 (Атлас)")

# Рис 4.2: Сравнение сигналов (До и После фильтрации)
plt.figure(figsize=(12, 8))
roi_idx = 0 # Индекс первой зоны для примера

plt.subplot(2, 1, 1)
plt.plot(raw_time_series[:, roi_idx], color='gray', alpha=0.6, label='Raw Signal')
plt.title(f'Raw BOLD Signal (ROI: {aal_atlas["labels"][roi_idx]})')
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(X[:, roi_idx], color='darkblue', label='Filtered & Standardized')
plt.title('Signal after Detrending, Bandpass Filtering (0.01-0.1 Hz) and Z-scoring')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'figure_4_2_filtering.png'), dpi=300)
print("✓ Сохранен Рис 4.2 (Фильтрация)")

# Рис 4.3: Тепловая карта (Heatmap)
plt.figure(figsize=(14, 7))
sns.heatmap(X.T, cmap='RdBu_r', center=0, cbar_kws={'label': 'Z-score'})
plt.title('Final Feature Matrix (116 Regions x Time)', fontsize=14)
plt.xlabel('Time (Volumes)')
plt.ylabel('Brain Regions (AAL Index)')
plt.savefig(os.path.join(output_dir, 'figure_4_3_heatmap.png'), dpi=300)
print("✓ Сохранен Рис 4.3 (Матрица признаков)")

plt.show()

print("\n" + "="*70)
print("Скрипт завершен успешно! Все данные и графики в папке /output")
print("="*70)