import tensorflow as tf
import pandas as pd
import uproot
import numpy as np
from matplotlib import pyplot as plt

#########################
# Config
training_file = "ml_inputs.root"
training_tree = "dummy_tree"
SHUFFLE_BUFFER = 500
BATCH_SIZE = 2

# Load data
print('Reading training data')
in_file = uproot.open(training_file)
in_tree = in_file.get(training_tree)
df_all = in_tree.arrays(library="pd")

df_all = df_all.drop(df_all[df_all.pT < 0].index).sample(frac=1)

# the test data will be 10% (0.1) of the entire data
test_size = int(len(df_all) * 0.1)

print(test_size)

train_df = df_all.iloc[:test_size, :].copy()
test_df = df_all.iloc[test_size:, :].copy()

print(train_df.head(5))

target = train_df.pop('target')

normalizer = tf.keras.layers.Normalization(axis=-1)
normalizer.adapt(train_df)

model = tf.keras.Sequential([
    normalizer,
    tf.keras.layers.Dense(10, activation='relu'),
    tf.keras.layers.Dense(10, activation='relu'),
    tf.keras.layers.Dense(1, activation='sigmoid')
])

model.compile(optimizer='adam',
              loss=tf.keras.losses.BinaryCrossentropy(),
              metrics=['accuracy'])

model.fit(train_df, target, epochs=15, batch_size=BATCH_SIZE)

test_sig = test_df.loc[test_df['target'] == 1]
test_bg = test_df.loc[test_df['target'] == 0]

tgs = test_sig.pop('target')
tgb = test_bg.pop('target')

print(test_sig.head(5))

test_target = test_df.pop('target')
test_loss, test_acc = model.evaluate(test_df, test_target)

pred_sig = list(model.predict(test_sig))
pred_bg = list(model.predict(test_bg))

probs_sig = pd.DataFrame(pred_sig)
probs_bg = pd.DataFrame(pred_bg)

# Let's look at the outputs of the linear and NN discriminants
fig, ax = plt.subplots(figsize=(5, 5))

probs_sig.plot(ax=ax, kind='hist', bins=100,
               title='predicted probabilities', alpha=0.7, label='sig')
probs_bg.plot(ax=ax, kind='hist', bins=100,
              title='predicted probabilities', alpha=0.7, label='bg')
ax.legend(frameon=False, prop={'size': 16})
fig.savefig('cippa.pdf')
